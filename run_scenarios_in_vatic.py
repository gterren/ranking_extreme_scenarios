import csv, pickle, sys, bz2, os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from mpi4py import MPI
from datetime import datetime

from vatic.data.loaders import GridLoader, load_input, RtsLoader, T7kLoader
from vatic.engines import Simulator

from typing import Tuple, Dict, Mapping, Optional, Union, Set, Iterable
from datetime import datetime
from pathlib import Path

# Define user paths
path_to_outputs   = r'/home/gterren/orfeus/outputs/Texas-7k/10_coupled'
path_to_scenarios = r'/home/gterren/orfeus/scenarios/clnSim'

# Define cohorts
cohort_ = ['Load', 'Wind', 'Solar']

# Specify grit to laod
_grid = T7kLoader()

# Simulation parameters
mipgap         = 0.02
reserve_factor = 0.10
sced_horizon   = 1
output_detail  = 2
verbosity      = 0
create_plots   = False
renew_costs    = None

# Define the number of scenarios to run
N_scen   = 1000
num_days = 1
coupled  = True

# Define date
day = sys.argv[1]
print(day)

# Crate directory if it does not exist
if not os.path.exists(path_to_outputs + "/{}".format(day)):
   os.makedirs(path_to_outputs + "/{}".format(day))

# Get MPI node information
def _get_node_info(verbose = True):
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    name = MPI.Get_processor_name()
    if verbose:
        print('>> MPI: Name: {} Rank: {} Size: {}'.format(name, rank, size) )
    return int(rank), int(size), comm

def _load_by_bus(self, start_date: Optional[datetime] = None,
                       end_date:   Optional[datetime] = None,
                       load_actls: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Parse forecast and actual load demands from zone to bus level."""
    load_fcsts = self.get_forecasts('Load', start_date, end_date)

    site_dfs = dict()
    for zone, zone_df in self.bus_df.groupby('Area'):
        area_total_load = zone_df['MW Load'].sum()

        for bus_name, bus_load in zip(zone_df['Bus Name'], zone_df['MW Load']):
            site_df            = pd.DataFrame({'fcst': load_fcsts[str(zone)],'actl': load_actls[str(zone)]}) * bus_load / area_total_load
            site_df.columns    = pd.MultiIndex.from_tuples([(x, bus_name) for x in site_df.columns])
            site_dfs[bus_name] = site_df

    return pd.concat(site_dfs.values(), axis = 1).sort_index(axis = 1)

def _create_scenario_timeseries(self, scen_dfs:   Mapping[str, pd.DataFrame],
                                      start_date: datetime,
                                      end_date:   datetime,
                                      scenario:   int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Get asset values corresponding to a scenario for actuals."""
    use_scens = {asset_type: scen_df.loc[scenario] for asset_type, scen_df in scen_dfs.items()}
    use_scens = {asset_type: scen_df.loc[scen_df.index >= start_date] for asset_type, scen_df in use_scens.items()}
    use_scens = {asset_type: scen_df.loc[scen_df.index <= end_date] for asset_type, scen_df in use_scens.items()}

    # consolidate scenario data for the renewable generators as the actuals
    gen_scens = pd.concat([use_scens['Wind'], use_scens['Solar']], axis = 1)

    gen_scens.columns = pd.MultiIndex.from_tuples([('actl', asset_name) for asset_name in gen_scens.columns])
    # consolidate forecasted output values for the renewable generators
    gen_df = pd.concat([self.get_forecasts(asset_type, start_date, end_date) for asset_type in self.timeseries_cohorts], axis = 1)

    # create one matrix with both forecasted and actual renewable outputs
    gen_df.columns = pd.MultiIndex.from_tuples([('fcst', asset_name) for asset_name in gen_df.columns])

    gen_df = pd.concat([gen_df, gen_scens], axis = 1)

    # for asset_type in self.no_scenario_renews:
    #     new_actuals = self.get_actuals(asset_type, start_date, end_date).resample('H').mean()
    #
    #     for asset_name, asset_actuals in new_actuals.iteritems():
    #         gen_df['actl', asset_name] = asset_actuals

    # TODO: investigate cases where forecasts are missing for assets
    #      with scenario values
    for gen in gen_df.actl.columns:
        if gen not in gen_df.fcst.columns:
            gen_df.drop(('actl', gen), axis = 1, inplace = True)

    assert (sorted(gen_df['fcst'].columns) == sorted(gen_df['actl'].columns)), ("Mismatching sets of assets with forecasts and with actuals!")

    demand_df = _load_by_bus(self, start_date, end_date, load_actls = use_scens['Load'])

    return gen_df.sort_index(axis = 1), demand_df

def _load_scenarios(self, scen_dir:    Union[str, Path],
                          dates:       Iterable[datetime],
                          scenarios:   Iterable[int],
                          asset_types: Tuple[str] = ('Load', 'Wind', 'Solar')) -> Dict[str, pd.DataFrame]:

    """Parse a set of scenarios saved to file."""
    scens = {asset_type: dict() for asset_type in asset_types}

    for scen_day in dates:
        out_file = Path(scen_dir, "clnSim_{}.p.gz".format(scen_day.strftime('%Y%m%d')))

        with bz2.BZ2File(out_file, 'r') as f:
            day_scens = pickle.load(f)

        for asset_type in asset_types:
            scens[asset_type][scen_day] = day_scens[asset_type].iloc[scenarios]

    scen_dfs = {asset_type: pd.concat(scens[asset_type].values(), axis = 1).stack() for asset_type in asset_types}

    for asset_type in asset_types:
        scen_dfs[asset_type].index = pd.MultiIndex.from_tuples([(scenario, t + self.utc_offset) for scenario, t in scen_dfs[asset_type].index], names = ('Scenario', 'Time'))

        if asset_type == 'Wind':
            scen_dfs[asset_type] = self.map_wind_generators(scen_dfs[asset_type])
        if asset_type == 'Solar':
            scen_dfs[asset_type] = self.map_solar_generators(scen_dfs[asset_type])

    return scen_dfs

# Process bus_detail.csv to get LMPs
def _get_lmps(buses_lmps_):
    buses_lmps_p_ = np.stack(buses_lmps_.index.to_numpy())
    buses_lmps_   = buses_lmps_.to_numpy()
    buses_hours_  = buses_lmps_p_[:, 1].astype(int)
    return np.stack([buses_lmps_[buses_hours_ == hour] for hour in range(24)])

# Process Vatic inputs to get actuals and forecast load, solar, and wind
def _fcst_and_actl(gen_, load_):
    gen_names_   = gen_.columns
    gen_values_  = gen_.to_numpy()[:24, :]
    load_names_  = load_.columns
    load_values_ = load_.to_numpy()[:24, :]
    # Find Wind and Solar generators
    idx_ = []
    for name in gen_names_:
        if name.find('Wind') > -1.:
            idx_.append(0)
        if name.find('Solar') > -1.:
            idx_.append(1)
    idx_ = np.array(idx_)
    # Aggregate assets in cohort
    load_values_  = np.sum(load_values_, axis = 1)
    solar_values_ = np.sum(gen_values_[:, idx_ == 1.], axis = 1)
    wind_values_  = np.sum(gen_values_[:, idx_ == 0.], axis = 1)
    net_values_   = load_values_ - (solar_values_ + wind_values_)
    return np.stack([load_values_, solar_values_, wind_values_, net_values_])

# Define starting and ending date
start_day = pd.Timestamp(day, tz = 'utc')
end_day   = start_day + pd.Timedelta(days = 1)

# MPI job variables
i_job, N_jobs, comm = _get_node_info()
print(i_job, N_jobs)
N_jobs_per_batches = int(N_scen/N_jobs)
print(N_jobs_per_batches)
idx_scen_ = [np.linspace(0, N_scen - 1, N_scen, dtype = int)[i_batch*N_jobs_per_batches:(i_batch + 1)*N_jobs_per_batches] for i_batch in range(N_jobs)]
print(len(idx_scen_))

simMap_    = pd.read_csv(path_to_scenarios + '/SimMap.csv')
simDate_   = np.array([int(date) for date in simMap_.columns])
simMap_    = simMap_.to_numpy()
simDate_p_ = pd.to_datetime(simDate_-70*365-18, unit='D')
index      = np.where(simDate_p_ == day)[0][0]
print(day, index, simDate_[index], simDate_p_[index], simMap_.shape)

for i_scen in idx_scen_[i_job]:

    # Define Scenario index for consecutive operational day
    idx_scen_1 = i_scen

    if coupled: idx_scen_2 = simMap_[i_scen, index]
    else:       idx_scen_2 = i_scen

    # Load Scenarios for each consecutive operational day
    scen_1_ = _load_scenarios(_grid, path_to_scenarios, [start_day], [idx_scen_1], cohort_)
    scen_2_ = _load_scenarios(_grid, path_to_scenarios, [end_day], [idx_scen_2], cohort_)

    # # Change the date and scenarios index within assets
    # for assets in cohort_:
    #     scen_2_index_ = scen_2_[assets].index.tolist()
    #     for i_hour in range(scen_2_[assets].index.values.shape[0]):
    #         scen_2_index_prime_    = list(scen_2_index_[i_hour])
    #         scen_2_index_prime_[0] = idx_scen_1
    #         scen_2_index_prime_[1] = scen_2_[assets].index.values[i_hour][1].replace(day = 2)
    #         scen_2_index_[i_hour]  = tuple(scen_2_index_prime_)
    #     scen_2_[assets].index = scen_2_index_
    # # Change the date and scenarios index within each asset
    # for assets in cohort_:
    #     for asset in scen_2_[assets].keys():
    #         scen_2_index_ = scen_2_[assets][asset].index.tolist()
    #         for i_hour in range(scen_2_[assets][asset].index.values.shape[0]):
    #             scen_2_index_prime_    = list(scen_2_index_[i_hour])
    #             scen_2_index_prime_[0] = idx_scen_1
    #             scen_2_index_prime_[1] = scen_2_[assets][asset].index.values[i_hour][1].replace(day = 2)
    #             scen_2_index_[i_hour]  = tuple(scen_2_index_prime_)
    #         scen_2_[assets][asset].index = scen_2_index_

    # Merge both Scenarios
    scen_ = {assets: pd.concat([scen_1_[assets], scen_2_[assets]], axis = 0) for assets in cohort_}

    gen_, load_ = _create_scenario_timeseries(_grid, scen_, start_day, end_day + pd.Timedelta(days = 1), idx_scen_1)
    print(gen_.shape, load_.shape)

    gen_p_  = gen_.drop(gen_.index[-24:])
    load_p_ = load_.drop(load_.index[-24:])

    # Simulation Variables
    start_date    = start_day.date()
    ruc_file      = None
    template_data = _grid.template
    gen_data      = gen_p_
    load_data     = load_p_

    # Define configuration of Vatic engine simulator
    _simulator = Simulator(template_data, gen_data, load_data, out_dir                        = None,
                                                               start_date                     = start_date,
                                                               num_days                       = 1,
                                                               solver                         = 'gurobi',
                                                               solver_options                 = {'Threads': 1},
                                                               run_lmps                       = True,
                                                               mipgap                         = mipgap,
                                                               reserve_factor                 = reserve_factor,
                                                               prescient_sced_forecasts       = False,
                                                               ruc_prescience_hour            = 0,
                                                               ruc_execution_hour             = 16,
                                                               ruc_every_hours                = 24,
                                                               ruc_horizon                    = 48,
                                                               sced_horizon                   = sced_horizon,
                                                               enforce_sced_shutdown_ramprate = False,
                                                               no_startup_shutdown_curves     = False,
                                                               output_detail                  = output_detail,
                                                               init_ruc_file                  = ruc_file,
                                                               verbosity                      = verbosity,
                                                               output_max_decimals            = 4,
                                                               create_plots                   = create_plots,
                                                               renew_costs                    = renew_costs,
                                                               save_to_csv                    = False)

    # Run simulation
    report_ = _simulator.simulate()

    # with open(path_to_outputs + "/{}".format(day) + r"/report_d{}_s{}.pkl".format(day, i_scen), 'wb') as _file:
    #     pickle.dump(report_, _file, protocol = pickle.HIGHEST_PROTOCOL)

    dict_keys(['hourly_summary', 'runtimes', 'total_runtime', 'ruc_summary', 'thermal_detail', 'renew_detail', 'daily_commits', 'bus_detail', 'line_detail'])
    summary_ = np.swapaxes(report_['hourly_summary'].to_numpy()[:, :-3], 0, 1)
    LMPs_    = _get_lmps(report_['bus_detail']['LMP'])
    actl_    = _fcst_and_actl(gen_p_['actl'], load_p_['actl'])
    fcst_    = _fcst_and_actl(gen_p_['fcst'], load_p_['fcst'])
    LMPs_p_  = np.swapaxes(np.concatenate((np.median(LMPs_, axis = 1)[:, np.newaxis],
                                           np.mean(LMPs_, axis = 1)[:, np.newaxis],
                                           np.min(LMPs_, axis = 1)[:, np.newaxis],
                                           np.max(LMPs_, axis = 1)[:, np.newaxis]), axis = 1), 0, 1)

    with open(path_to_outputs + "/{}".format(day) + r"/d{}_s{}.pkl".format(day, i_scen), 'wb') as _file:
        pickle.dump([summary_, LMPs_p_, actl_, fcst_], _file, protocol = pickle.HIGHEST_PROTOCOL)


