import numpy
import pandas
from natcap.invest import datastack

# For each subregion
#    For each timestep
#        For each age/stage

def get_subregion_names(dataframe):
    return


def lobster():
    paramset = datastack.extract_parameter_set(
        '../../invest/data/invest-sample-data/spiny_lobster_belize.invs.json')
    args = paramset.args.copy()

    alpha = numpy.float64(args['alpha'])
    beta = numpy.float64(args['beta'])
    n_sexes = numpy.int64(args['sexsp'])
    n_init_recruits = numpy.int64(args['total_init_recruits'])
    n_timesteps = numpy.int64(args['total_timesteps'])
    unit_price = numpy.float64(args['unit_price'])

    # parsed population parameters
    population_params = pandas.read_csv(param_csv_path)
    per_subregion_params = {}
    per_stage_params = {}
    for subregion in population_params.columns:
        if subregion.lower() == 'age_area':
            continue

        per_subregion_params[subregion] = {}

        if subregion.startwith('Unnamed:'):
            break

        for row in population_params.iterrows():
            if row.isna():
                break

            stage_name = row[Age_Area]

            stage_parameters = {}
            stage_parameters['Mortality'] = row[subregion]
            stage_parameters['Maturity'] = row.
            stage_parameters['ExploitationFraction'] = row.ExploitationFraction
            stage_parameters['LarvalDispersal'] = row.LarvalDispersal

            per_subregion_params[subregion][stage_name] = stage_parameters

    parameter_start_index = float('inf')
    for row_index, row in population_params.iterrows():
        if row.isna():
            parameter_start_index = row_index
        else:
            if row_index < parameter_start_index:
                continue

        parameter_name = row['Age_Area']

        for subregion_name in population_params.columns:
            if subregion_name.lower() == 'age_area':
                continue

            if numpy.isnan(row[subregion_name]):
                continue

            if subregion_name not in per_subregion_params:
                per_subregion_params[subregion_name] = {}

            per_subregion_params[subregion_name][parameter_name] = row[subregion_name]

    total_spawners = [0]  # indexed by timestep
    total_recruits = [n_init_recruits]  # indexed by timestep

    for timestep in range(1, n_timesteps + 1):
        pass




