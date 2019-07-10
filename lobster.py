import pprint

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

    pprint.pprint(args)

    alpha = numpy.float64(args['alpha'])
    beta = numpy.float64(args['beta'])
    if args['sexsp'].lower() == 'yes':
        n_sexes = 2
    else:
        n_sexes = 1
    n_init_recruits = numpy.int64(args['total_init_recruits'])
    n_timesteps = numpy.int64(args['total_timesteps'])
    unit_price = numpy.float64(args['unit_price'])

    # parsed population parameters
    population_params = pandas.read_csv(args['population_csv_path'])
    per_subregion_params = {}
    per_stage_params = {}
    for subregion in population_params.columns:
        if subregion.lower() == 'age_area':
            continue

        if subregion.startswith('Unnamed:'):
            break

        per_subregion_params[subregion] = {'stages': {}}

        for row_index, row in population_params.iterrows():
            if all(row.isna()):
                break

            stage_name = row['Age_Area']

            stage_parameters = {}
            stage_parameters['Mortality'] = row[subregion]
            stage_parameters['Maturity'] = row.Maturity
            stage_parameters['VulnFishing'] = row.VulnFishing
            stage_parameters['Weight'] = row.Weight

            per_subregion_params[subregion]['stages'][stage_name] = stage_parameters

    parameter_start_index = float('inf')
    for row_index, row in population_params.iterrows():
        if all(row.isna()):
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

    pprint.pprint(per_subregion_params)

    total_spawners = [0]  # indexed by timestep
    total_recruits = [n_init_recruits]  # indexed by timestep

    for timestep in range(1, n_timesteps + 1):
        pass


if __name__ == '__main__':
    lobster()
