import collections
import pprint

import numpy
import pandas
from natcap.invest import datastack


def beverton_holt(alpha, beta, spawners):
    return (alpha*spawners) / (1 + (beta*spawners))


def lobster():
    paramset = datastack.extract_parameter_set(
        '../../invest/data/invest-sample-data/spiny_lobster_belize.invs.json')
    args = paramset.args.copy()

    pprint.pprint(args)

    # lobster is age-based.
    #alpha = numpy.float64(args['alpha'])
    #beta = numpy.float64(args['beta'])
    # n_init_recruits = numpy.int64(args['total_init_recruits'])
    n_timesteps = numpy.int64(args['total_timesteps'])
    unit_price = numpy.float64(args['unit_price'])

    n_init_recruits = 4686959
    alpha = numpy.float64(1000)
    beta = numpy.float64(0.00000016069)

    if args['sexsp'].lower() == 'yes':
        n_sexes = 2.0
    else:
        n_sexes = 1.0

    # parsed population parameters
    population_params = pandas.read_csv(args['population_csv_path'])
    per_subregion_params = {}
    per_stage_params = {}
    for subregion in population_params.columns:
        if subregion.lower() == 'age_area':
            continue

        if subregion.startswith('Unnamed:'):
            break

        per_subregion_params[subregion] = {'stages': collections.OrderedDict()}

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
    n_stages = len(per_subregion_params.values()[0]['stages'])

    # survival[subregion][stage_index]
    survival = collections.defaultdict(list)
    for subregion, subregion_params in per_subregion_params.items():
        exploitation = subregion_params['ExploitationFraction']

        for stage_index, stage in enumerate(subregion_params['stages'].values()):
            survival[subregion].append(
                stage['Mortality']*(1-(exploitation * stage['VulnFishing'])))

    # populations[subregion][timestep][stage_index]
    populations = collections.defaultdict(lambda: collections.defaultdict(list))
    for timestep in range(0, n_timesteps + 1):
        total_recruits = 0
        total_population = 0

        # Set initial conditions (different from normal modelling)
        if timestep == 0:
            for subregion, subregion_params in per_subregion_params.items():
                larval_dispersal = subregion_params['LarvalDispersal']
                for stage_index, stage in enumerate(subregion_params['stages'].values()):
                    if stage_index == 0:
                        population = (n_init_recruits/float(n_sexes)) * larval_dispersal
                    elif stage_index < (n_stages -1):
                        population = populations[subregion][0][stage_index-1] * survival[subregion][stage_index-1]
                    else:
                        population = (
                            (populations[subregion][0][stage_index-1] * survival[subregion][stage_index-1]) /
                            (1 - survival[subregion][stage_index]))

                    populations[subregion][timestep].append(population)
        else:
            for subregion, subregion_params in per_subregion_params.items():
                larval_dispersal = subregion_params['LarvalDispersal']
                migration = 1  # No migration for now.

                for stage_index, stage in enumerate(subregion_params['stages'].values()):
                    if stage_index == 0:
                        population = beverton_holt(alpha, beta, spawners=n_init_recruits/n_sexes)
                    elif stage_index < (n_stages - 1):
                        population = populations[subregion][timestep-1][stage_index-1] * survival[subregion][stage_index-1]
                    else:  # We're at max stage
                        population = (
                            (populations[subregion][timestep-1][stage_index-1] * survival[subregion][stage_index-1]) +
                            (populations[subregion][timestep-1][stage_index] * survival[subregion][stage_index]))

                    # Can't have fewer than 0 individuals
                    # Can't have fractional individuals.
                    population = max(0, population)

                    populations[subregion][timestep].append(population)

    # produce 1 table per subregion












if __name__ == '__main__':
    lobster()
