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
    model_type = 'age'
    #alpha = numpy.float64(args['alpha'])
    #beta = numpy.float64(args['beta'])
    # n_init_recruits = numpy.int64(args['total_init_recruits'])

    # These are the parameters that Jess and Jodie want to use, but their
    # migrations spreadsheet isn't referencing their modified sheet.
    n_init_recruits = 4686959
    alpha = numpy.float64(1000)
    beta = numpy.float64(0.00000016069)
    recruitment = beverton_holt_1

    # These are the values in the Model_Lobster sheet
    #n_init_recruits = 4686959
    #alpha = numpy.float64(5770000)
    #beta = numpy.float64(2885000)
    #recruitment = beverton_holt_2


    model(args)

def shrimp():
    # Uses fixed recruitment
    # Stage-based model
    # No migration
    paramset = datastack.extract_parameter_set(
        '../../invest/data/invest-sample-data/white_shrimp_galveston_bay.invs.json')
    args = paramset.args.copy()

    model(args)


def model(args):
    n_timesteps = numpy.int64(args['total_timesteps'])
    n_init_recruits = numpy.int64(args['total_init_recruits'])

    if args['sexsp'].lower() == 'yes':
        n_sexes = 2.0
    else:
        n_sexes = 1.0

    if args['population_type'].lower().startswith('stage'):
        model_type = 'stage'
    else:
        model_type = 'age'

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

            try:
                # This the stage name header
                stage_name = row['Age_Area']
            except KeyError:
                # shrimp uses this header instead.
                stage_name = row['Stage']

            stage_parameters = {}
            stage_parameters['Mortality'] = row[subregion]

            if args['recruitment_type'].lower() in ('beverton-holt', 'ricker'):
                stage_parameters['Maturity'] = row.Maturity
            else:
                stage_parameters['Maturity'] = 1

            stage_parameters['VulnFishing'] = row.VulnFishing
            stage_parameters['Weight'] = row.Weight
            stage_parameters['Duration'] = row.Duration

            per_subregion_params[subregion]['stages'][stage_name] = stage_parameters

    parameter_start_index = float('inf')
    for row_index, row in population_params.iterrows():
        if all(row.isna()):
            parameter_start_index = row_index
        else:
            if row_index < parameter_start_index:
                continue

        try:
            # This the parameter name header
            parameter_name = row['Age_Area']
        except KeyError:
            # shrimp uses this header instead.
            parameter_name = row['Stage']

        for subregion_name in population_params.columns:
            if subregion_name.lower() in ('age_area', 'stage'):
                continue

            if numpy.isnan(row[subregion_name]):
                continue

            if subregion_name not in per_subregion_params:
                per_subregion_params[subregion_name] = {}

            per_subregion_params[subregion_name][parameter_name] = row[subregion_name]

    total_spawners = [0]  # indexed by timestep
    total_recruits = [n_init_recruits]  # indexed by timestep
    n_stages = len(per_subregion_params.values()[0]['stages'])
    subregions = list(per_subregion_params.keys())

    # survival[subregion][stage_index]
    survival = collections.defaultdict(list)
    for subregion, subregion_params in per_subregion_params.items():
        try:
            exploitation = subregion_params['ExploitationFraction']
        except KeyError:
            import pdb; pdb.set_trace()

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

    def population_after_migration(this_subregion, timestep, stage_index):
        return sum(
            (populations[other_subregion][timestep][stage_index] *
             migration[stage_name][other_subregion][this_subregion])
            for other_subregion in subregions)

    # populations[subregion][timestep][stage_index]
    populations = collections.defaultdict(lambda: collections.defaultdict(list))
    for timestep in range(0, n_timesteps + 1):
        spawners = 0

        for subregion, subregion_params in per_subregion_params.items():
            larval_dispersal_per_sex = subregion_params['LarvalDispersal'] / n_sexes

            for stage_index, (stage_name, stage) in enumerate(subregion_params['stages'].items()):
                assert stage_index == int(stage_name)  # should be true for test data.


                if model_type == 'age':
                    if stage_index == 0:
                        population = total_recruits[timestep] * larval_dispersal_per_sex

                    elif stage_index < (n_stages -1):
                        # disallow migration in the initial timestep.
                        if stage_name in migration and timestep > 0:
                            population = population_after_migration(subregion, timestep-1, stage_index-1)

                        else:
                            # Don't allow the population index to go below 0.
                            # If we're setting the initial condition, this will be:
                            #    populations[subregion][0][stage_index-1]
                            #
                            # Otherwise, this will be
                            #    populations[subregion][timestep-1][stage_index-1]
                            population = populations[subregion][max(0, timestep-1)][stage_index-1]
                        population *= survival[subregion][stage_index-1]

                    else:
                        survival_from_previous_stage = survival[subregion][stage_index-1]
                        survival_from_final_stage = survival[subregion][stage_index]

                        # Initial conditions are a bit special for the final stage of a subregion.
                        # Migration is disallowed for this timestep.
                        if timestep == 0:
                            population = (
                                (populations[subregion][0][stage_index-1] * survival_from_previous_stage) /
                                (1 - survival_from_final_stage))

                        # Migration can only happen if it's allowed for this stage
                        # and we're not at timestep 0.
                        if stage_name in migration and timestep > 0:
                            prev_stage_population = population_after_migration(
                                subregion, timestep-1, stage_index-1) * survival_from_previous_stage

                            final_stage_population = population_after_migration(
                                subregion, timestep-1, stage_index)* survival_from_final_stage

                            population = prev_stage_population + final_stage_population

                        # No migration, we're not at timestep 0 and we're at the final stage.
                        else:
                            population = (
                                (populations[subregion][timestep-1][stage_index-1] * survival_from_previous_stage) +
                                (populations[subregion][timestep-1][stage_index] * survival_from_final_stage))

                # Stage-base populations have a completely different model structure.
                elif model_type == 'stage':
                    if timestep == 0:
                        if stage_index == 0:
                            population = total_recruits[timestep] * larval_dispersal_per_sex
                        else:
                            population = 1
                    else:
                        # P_asx in the User's Guide.  The probability of
                        # surviving from natural and fishing mortality and
                        # *staying* in the same stage for this sex and subregion.
                        stage_survival = survival[subregion][stage_index]
                        probability_of_staying_in_same_stage = (
                            stage_survival * (
                                (1-stage_survival**(stage['Duration']-1)) /
                                (1-stage_survival**(stage['Duration']))))

                        if stage_index == 0:
                            population = ((
                                population_after_migration(
                                    subregion, timestep-1, stage_index) *
                                probability_of_staying_in_same_stage) + (
                                    total_recruits[timestep]))
                        else:
                            # G_asx in the User's Guide.  The probability of
                            # surviving from natural and fishing mortality and
                            # *growing* into the next stage for each sex and
                            # subregion.
                            probability_of_growing_to_next_stage = (
                                (stage_survival**stage['Duration'] * (1-stage_survival)) /
                                (1-stage_survival**stage['Duration']))
                            population = (
                                (population_after_migration(
                                    subregion, timestep-1, stage_index-1) *
                                 probability_of_growing_to_next_stage) +
                                (population_after_migration(
                                    subregion, timestep-1, stage_index) *
                                 probability_of_staying_in_same_stage))
                else:
                    raise AssertionError('Invalid model type')

                # Can't have fewer than 0 individuals
                population = max(0, population)

                populations[subregion][timestep].append(population)
                spawners += stage['Maturity'] * stage['Weight'] * population

        total_spawners.append(spawners)
        total_recruits.append(max(0, recruitment(alpha, beta, spawners)))

    # For now, just produce a dataframe for a single subregion for all timesteps.
    out_df = pandas.DataFrame.from_dict(populations['1'], orient='index')
    print(out_df)


if __name__ == '__main__':
    #lobster()
    shrimp()
