import math
import os
import glob
import collections
import pprint

import numpy
import pandas
from natcap.invest import datastack


def beverton_holt_1(alpha, beta, spawners):
    return float(alpha*spawners) / (1. + (beta*spawners))


def beverton_holt_2(alpha, beta, spawners):
    return float(alpha * spawners) / (beta + spawners)


def ricker(alpha, beta, spawners):
    return (alpha*spawners*math.e**(-beta*spawners))


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

    # These are the parameters that Jess and Jodie want to use, but their
    # migrations spreadsheet isn't referencing their modified sheet.
    n_init_recruits = 4686959
    alpha = numpy.float64(1000)
    beta = numpy.float64(0.00000016069)
    recruitment = beverton_holt_1

    # These are the values in the Model_Lobster sheet
    n_init_recruits = 4686959
    alpha = numpy.float64(5770000)
    beta = numpy.float64(2885000)
    recruitment = beverton_holt_2

    if args['sexsp'].lower() == 'yes':
        n_sexes = 2.0
    else:
        n_sexes = 1.0

    # parsed population parameters
    population_params = pandas.read_csv(args['population_csv_path'])
    per_subregion_params = collections.OrderedDict()
    per_stage_params = collections.OrderedDict()
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

    # migration[stage_name][from_subregion][to_subregion] = float
    migration = {}
    for csv_filename in glob.glob(os.path.join(args['migration_dir'], '*.csv')):
        stage_name = os.path.splitext(os.path.basename(csv_filename))[0].split('_')[-1]
        temp_dict = pandas.DataFrame.to_dict(
            pandas.read_csv(csv_filename, index_col=0), orient='dict')

        # Pandas' to_dict() reads in the sink subregions as integers if they
        # are integers, so we need to make sure that they are all strings so we
        # can index the dict properly.
        cast_keys_dict = {}
        for source_subregion, sink_dict in temp_dict.items():
            nested_dict = {}
            for sink_subregion, migration_proportion in sink_dict.items():
                nested_dict[str(sink_subregion)] = migration_proportion

            cast_keys_dict[source_subregion] = nested_dict
        migration[stage_name] = cast_keys_dict

    # For keeping track of spawners and recruits
    # Indexed by timestep.
    total_spawners = [0]
    total_recruits = [n_init_recruits]

    # populations[subregion][timestep][stage_index]
    populations = collections.defaultdict(lambda: collections.defaultdict(list))
    for timestep in range(0, n_timesteps + 1):
        spawners = 0

        # Set initial conditions (different from normal modelling)
        if timestep == 0:
            for subregion, subregion_params in per_subregion_params.items():
                larval_dispersal_per_sex = subregion_params['LarvalDispersal'] / n_sexes

                for stage_index, stage in enumerate(subregion_params['stages'].values()):
                    if stage_index == 0:
                        population = n_init_recruits * larval_dispersal_per_sex
                    elif stage_index < (n_stages -1):
                        population = populations[subregion][0][stage_index-1] * survival[subregion][stage_index-1]
                    else:
                        population = (
                            (populations[subregion][0][stage_index-1] * survival[subregion][stage_index-1]) /
                            (1 - survival[subregion][stage_index]))

                    populations[subregion][timestep].append(population)
                    spawners += stage['Maturity'] * stage['Weight'] * population
        else:
            for subregion, subregion_params in per_subregion_params.items():
                larval_dispersal_per_sex = subregion_params['LarvalDispersal'] / n_sexes

                for stage_index, (stage_name, stage) in enumerate(subregion_params['stages'].items()):
                    assert stage_index == int(stage_name)  # should be true for test data.

                    if stage_index == 0:
                        population = total_recruits[timestep] * larval_dispersal_per_sex
                    elif stage_index < (n_stages - 1):

                        # If there is defined migration for this stage, add that calculation in here.
                        if stage_name in migration:
                            population = 0
                            for other_subregion in per_subregion_params.keys():
                                prev_stage_population = populations[other_subregion][timestep-1][stage_index-1]
                                migration_proportion = migration[stage_name][other_subregion][subregion]
                                population += (prev_stage_population * migration_proportion)
                            population *= survival[subregion][stage_index-1]
                        else:
                            population = populations[subregion][timestep-1][stage_index-1] * survival[subregion][stage_index-1]

                    else:  # We're at max stage
                        if stage_name in migration:
                            prev_stage_population = 0
                            for other_subregion in per_subregion_params.keys():
                                prev_stage_population = populations[other_subregion][timestep-1][stage_index-1]
                                migration_proportion = migration[stage_name][other_subregion][subregion]
                                prev_stage_population += (prev_stage_population * migration_proportion)
                            prev_stage_population *= survival[subregion][stage_index-1]

                            final_stage_population = 0
                            for other_subregion in per_subregion_params.keys():
                                final_stage_population += (
                                    populations[other_subregion][timestep-1][stage_index] *
                                    migration[stage_name][other_subregion][subregion])
                            final_stage_population *= survival[subregion][stage_index]

                            population = prev_stage_population + final_stage_population
                        else:
                            population = (
                                (populations[subregion][timestep-1][stage_index-1] * survival[subregion][stage_index-1]) +
                                (populations[subregion][timestep-1][stage_index] * survival[subregion][stage_index]))

                    # Can't have fewer than 0 individuals
                    # Can't have fractional individuals.
                    population = max(0, population)

                    populations[subregion][timestep].append(population)
                    spawners += stage['Maturity'] * stage['Weight'] * population

        total_spawners.append(spawners)
        total_recruits.append(max(0, recruitment(alpha, beta, spawners)))

    # For now, just produce a dataframe for a single subregion for all timesteps.
    out_df = pandas.DataFrame.from_dict(populations['1'], orient='index')
    print(out_df)


if __name__ == '__main__':
    lobster()
