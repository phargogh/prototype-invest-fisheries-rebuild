import math
import os
import glob
import collections
import pprint
import logging

import numpy
import pandas
from natcap.invest import datastack


logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger(__name__)


def beverton_holt_1(args):
    alpha = args['alpha']
    beta = args['beta']
    spawners = args['spawners']
    return float(alpha*spawners) / (1. + (beta*spawners))

def beverton_holt_2(args):
    alpha = args['alpha']
    beta = args['beta']
    spawners = args['spawners']
    return float(alpha * spawners) / (beta + spawners)

def fecundity(args):
    per_subregion_params = args['per_subregion_params']
    populations = args['populations']
    recruits = 0
    for subregion, subregion_params in per_subregion_params.items():
        for sex_index in sexes:
            for stage_index, (stage_name, stage) in enumerate(subregion_params['stages'][sex_index].items()):
                recruits += (
                    populations[sex_index][subregion][timestep-1][stage_index-1] *
                    stage['Maturity'] * stage['Fecundity'])

    return recruits


def ricker(args):
    alpha = args['alpha']
    beta = args['beta']
    spawners = args['spawners']
    return (alpha*spawners*math.e**(-beta*spawners))

def fixed(args):
    return args['n_recruits']

def check_expected(actual, expected):
    assert round(expected) == expected, '%s != %s' % (actual, expected)


def lobster():
    paramset = datastack.extract_parameter_set(
        '../../invest/data/invest-sample-data/spiny_lobster_belize.invs.json')
    args = paramset.args.copy()
    args['total_init_recruits'] = 4686959  # to match spreadsheet

    LOGGER.info('Spiny Lobster - Sample Data')
    spawners, harvest = model(args, recruitment=beverton_holt_2)

    check_expected(spawners, 2847870)
    check_expected(harvest, 963451)


def lobster_jess():
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
    args['alpha'] = alpha
    args['beta'] = beta
    args['total_init_recruits'] = n_init_recruits

    # These are the values in the Model_Lobster sheet
    #n_init_recruits = 4686959
    #alpha = numpy.float64(5770000)
    #beta = numpy.float64(2885000)
    #recruitment = beverton_holt_2

    LOGGER.info('Spiny Lobster - Jess')
    model(args, recruitment=recruitment)

def shrimp():
    # Uses fixed recruitment
    # Stage-based model
    # No migration
    LOGGER.info('White Shrimp')
    paramset = datastack.extract_parameter_set(
        '../../invest/data/invest-sample-data/white_shrimp_galveston_bay.invs.json')
    args = paramset.args.copy()

    spawners, harvest = model(args, recruitment=fixed)

    check_expected(spawners, 2.16e11)  # fixed recruitment
    check_expected(harvest, 3096000)


def blue_crab():
    # Ricker recruitment
    # spawners are individuals
    # age-based population model
    # non-sex-specific population
    # harvest by individuals
    LOGGER.info('Blue Crab')
    paramset = datastack.extract_parameter_set(
        '../../invest/data/invest-sample-data/blue_crab_galveston_bay.invs.json')
    args = paramset.args.copy()

    spawners, harvest = model(args, recruitment=ricker)
    check_expected(spawners, 42644460)
    check_expected(harvest, 24798419)


def dungeness_crab():
    # Ricker recruitment
    # spawners are individuals
    # age-based population model
    # no migration
    # Sex-specific population
    LOGGER.info('Dungeness Crab')
    paramset = datastack.extract_parameter_set(
        '../../invest/data/invest-sample-data/dungeness_crab_hood_canal.invs.json')
    args = paramset.args.copy()

    spawners, harvest = model(args, recruitment=ricker)
    check_expected(spawners, 4051538)
    check_expected(harvest, 526987)


def model(args, recruitment):
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

    try:
        alpha = args['alpha']
    except KeyError:
        alpha = None

    try:
        beta = args['beta']
    except KeyError:
        beta = None

    try:
        n_recruits = args['total_recur_recruits']
    except KeyError:
        n_recruits = None


    # parsed population parameters
    population_params = pandas.read_csv(args['population_csv_path'])
    per_subregion_params = collections.OrderedDict()
    per_stage_params = collections.OrderedDict()
    sexes = list(range(0, int(n_sexes)))
    assert(len(sexes) == n_sexes)

    # Subregions may only be in columns 1+
    for subregion in population_params.columns[1:]:
        if subregion.startswith('Unnamed:'):
            break

        per_subregion_params[subregion] = {'stages': [collections.OrderedDict() for sex in sexes]}

        encountered_stage_names = set()
        for row_index, row in population_params.iterrows():
            if all(row.isna()):
                break

            # The stage name is always in the first column
            stage_name = row[0]
            if stage_name not in encountered_stage_names:
                _sex_index = 0
            else:
                if n_sexes == 1:
                    raise ValueError(
                        'Two stages have the same name but this '
                        'model is not sex-specific')

                # We're in a new sex of an existing stage.
                _sex_index = 1
            encountered_stage_names.add(stage_name)


            stage_parameters = {}
            stage_parameters['Mortality'] = row[subregion]

            if args['recruitment_type'].lower() in ('beverton-holt', 'ricker'):
                stage_parameters['Maturity'] = row.Maturity
            else:
                stage_parameters['Maturity'] = 1

            stage_parameters['VulnFishing'] = row.VulnFishing
            try:
                stage_parameters['Weight'] = row.Weight
            except AttributeError:
                if (args['harvest_units'].lower() == 'weight' or
                        args['spawn_units'].lower() == 'weight'):
                    raise ValueError('Weight vector required when harvest units are weight')
                stage_parameters['Weight'] = 1

            try:
                stage_parameters['Duration'] = row.Duration
            except AttributeError:
                # Duration column required in stage-based models.
                if model_type == 'stage':
                    raise

            per_subregion_params[subregion]['stages'][_sex_index][stage_name] = stage_parameters

    parameter_start_index = float('inf')
    for row_index, row in population_params.iterrows():
        if all(row.isna()):
            parameter_start_index = row_index
            continue
        else:
            if row_index < parameter_start_index:
                continue

        # Parameter name is always in the first column
        parameter_name = row[0]

        for subregion_name in population_params.columns[1:]:
            if numpy.isnan(row[subregion_name]):
                continue

            per_subregion_params[subregion_name][parameter_name] = row[subregion_name]

    stages = list(per_subregion_params.values()[0]['stages'][0].keys())
    subregions = list(per_subregion_params.keys())

    # survival[sex_index][subregion][stage_index]
    survival = [collections.defaultdict(list) for sex in sexes]
    for subregion, subregion_params in per_subregion_params.items():
        exploitation = subregion_params['ExploitationFraction']

        for sex_index in sexes:
            for stage_index, stage in enumerate(subregion_params['stages'][sex_index].values()):
                survival[sex_index][subregion].append(
                    stage['Mortality']*(1-(exploitation * stage['VulnFishing'])))

    # migration[stage_name][from_subregion][to_subregion] = float
    migration = {}
    if 'migration_dir' in args and args['migration_dir'] not in (None, ''):
        for csv_filename in glob.glob(os.path.join(args['migration_dir'], '*.csv')):
            stage_name = os.path.splitext(os.path.basename(csv_filename))[0].split('_')[-1]
            temp_df = pandas.read_csv(csv_filename, index_col=0)
            temp_dict = pandas.DataFrame.to_dict(
                temp_df, orient='dict')

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
    else:
        LOGGER.info('No migration provided')

    # For keeping track of spawners and recruits
    # Indexed by timestep.
    total_spawners = [0]
    total_recruits = [n_init_recruits]

    def population_after_migration(this_subregion, this_stage_name, timestep, stage_index, sex_index):
        return sum(
            (populations[sex_index][other_subregion][timestep][stage_index] *
             migration[this_stage_name][other_subregion][this_subregion])
            for other_subregion in subregions)

    # populations[sex_index][subregion][timestep][stage_index]
    populations = [collections.defaultdict(lambda: collections.defaultdict(list)) for sex in sexes]
    for timestep in range(0, n_timesteps + 1):
        spawners = 0

        for sex_index in sexes:
            for subregion, subregion_params in per_subregion_params.items():
                try:
                    larval_dispersal = subregion_params['LarvalDispersal']
                except KeyError:
                    # UG specifies that larvae disperse equally across all
                    # subregions unless the user specifies a different proportion.
                    larval_dispersal = 1. / len(subregions)
                larval_dispersal_per_sex = larval_dispersal / n_sexes

                for stage_index, (stage_name, stage) in enumerate(subregion_params['stages'][sex_index].items()):

                    if model_type == 'age':
                        if stage_index == 0:
                            population = total_recruits[timestep] * larval_dispersal_per_sex

                        elif stage_index < (len(stages) - 1):
                            # disallow migration in the initial timestep.
                            if stage_name in migration and timestep > 0:
                                population = population_after_migration(subregion, stage_name, timestep-1, stage_index-1, sex_index)

                            else:
                                # Don't allow the population index to go below 0.
                                # If we're setting the initial condition, this will be:
                                #    populations[sex_index][subregion][0][stage_index-1]
                                #
                                # Otherwise, this will be
                                #    populations[sex_index][subregion][timestep-1][stage_index-1]
                                population = populations[sex_index][subregion][max(0, timestep-1)][stage_index-1]
                            population *= survival[sex_index][subregion][stage_index-1]

                        else:
                            survival_from_previous_stage = survival[sex_index][subregion][stage_index-1]
                            survival_from_final_stage = survival[sex_index][subregion][stage_index]

                            # Initial conditions are a bit special for the final stage of a subregion.
                            # Migration is disallowed for this timestep.
                            if timestep == 0:
                                population = (
                                    (populations[sex_index][subregion][0][stage_index-1] * survival_from_previous_stage) /
                                    (1 - survival_from_final_stage))

                            # Migration can only happen if it's allowed for this stage
                            # and we're not at timestep 0.
                            elif stage_name in migration and timestep > 0:
                                prev_stage_population = population_after_migration(
                                    subregion, stage_name, timestep-1, stage_index-1, sex_index) * survival_from_previous_stageA

                                final_stage_population = population_after_migration(
                                    subregion, stage_name, timestep-1, stage_index, sex_index)* survival_from_final_stage

                                population = prev_stage_population + final_stage_population

                            # No migration, we're not at timestep 0 and we're at the final stage.
                            else:
                                population = (
                                    (populations[sex_index][subregion][timestep-1][stage_index-1] * survival_from_previous_stage) +
                                    (populations[sex_index][subregion][timestep-1][stage_index] * survival_from_final_stage))

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
                            stage_survival = survival[sex_index][subregion][stage_index]
                            probability_of_staying_in_same_stage = (
                                stage_survival * (
                                    (1-stage_survival**(stage['Duration']-1)) /
                                    (1-stage_survival**(stage['Duration']))))

                            if stage_index == 0:
                                if stage_name in migration:
                                    population = ((
                                        population_after_migration(
                                            subregion, stage_name, timestep-1, stage_index, sex_index) *
                                        probability_of_staying_in_same_stage) + (
                                            total_recruits[timestep]))
                                else:
                                    population = (
                                        populations[sex_index][subregion][timestep-1][stage_index] *
                                        probability_of_staying_in_same_stage +
                                        total_recruits[timestep])
                            else:
                                # G_asx in the User's Guide.  The probability of
                                # surviving from natural and fishing mortality and
                                # *growing* into the next stage for each sex and
                                # subregion.
                                probability_of_growing_to_next_stage = (
                                    (stage_survival**stage['Duration'] * (1-stage_survival)) /
                                    (1-stage_survival**stage['Duration']))

                                if stage_name in migration:
                                    population = (
                                        (population_after_migration(
                                            subregion, stage_name, timestep-1, stage_index-1, sex_index) *
                                         probability_of_growing_to_next_stage) +
                                        (population_after_migration(
                                            subregion, stage_name, timestep-1, stage_index, sex_index) *
                                         probability_of_staying_in_same_stage))
                                else:
                                    population = (
                                        (populations[sex_index][subregion][timestep-1][stage_index-1] *
                                         probability_of_growing_to_next_stage) +
                                        (populations[sex_index][subregion][timestep-1][stage_index] *
                                         probability_of_staying_in_same_stage))
                    else:
                        raise AssertionError('Invalid model type')

                    # Can't have fewer than 0 individuals
                    population = max(0, population)

                    populations[sex_index][subregion][timestep].append(population)

                    # This is where we handle spawners by weight (if requested)
                    stage_spawners = population * stage['Maturity']
                    if args['spawn_units'].lower() == 'weight':
                        stage_spawners*= stage['Weight']
                    spawners += stage_spawners

        total_spawners.append(spawners)
        recruitment_args = {
            'alpha': alpha,
            'beta': beta,
            'spawners': spawners,
            'n_recruits': n_recruits,
            'per_subregion_params': per_subregion_params,
            'populations': populations,
        }
        total_recruits.append(max(0, recruitment(recruitment_args)))

    # For now, just produce a dataframe for a single subregion for all timesteps.
    # This dataframe is for populations from subregion 1 (sex 0) only.
    out_df = pandas.DataFrame.from_dict(populations[0]['1'], orient='index')
    out_df.columns = stages  # reset stage names from indexes
    LOGGER.info('POPULATION SUBREGION 1 SEX 0')
    print(out_df)

    # This harvest dataframe is the harvest calculated for all subregions.
    # This would be an intermediate output.
    # harvest[subregion][timestep]
    harvest = {}
    for subregion in subregions:
        # exploitation is the 'harvest' row in the spreadsheet,
        # 'ExploitationFraction' in the population parameters table.
        exploitation = per_subregion_params[subregion]['ExploitationFraction']
        harvest[subregion] = []
        for timestep in range(0, n_timesteps+1):
            harvest_population = 0
            for sex_index in sexes:
                for stage_index, stage_name in enumerate(stages):
                    # vulnerability is the 'vulnerability' row in the spreadsheet,
                    # represented in populations input as VulnFishing
                    vulnerability = per_subregion_params[subregion]['stages'][sex_index][stage_name]['VulnFishing']
                    weight = per_subregion_params[subregion]['stages'][sex_index][stage_name]['Weight']

                    harvest_population += (
                        populations[sex_index][subregion][timestep][stage_index] *
                        exploitation *
                        vulnerability *
                        weight)

            harvest[subregion].append(harvest_population)
    harvest_df = pandas.DataFrame.from_dict(harvest, orient='columns')
    LOGGER.info('HARVEST')
    print(harvest_df)

    # This harvest value dataframe is from the final harvest numbers.
    if args['val_cont']:
        harvest_value = []
        total_harvest = 0
        for subregion in subregions:
            harvest_value.append({
                'subregion': subregion,
                'harvest': harvest[subregion][-1],
                'value': harvest[subregion][-1] * args['unit_price'] * args['frac_post_process'],
            })
            total_harvest += harvest[subregion][-1]

        harvest_value.append({
            'subregion': 'TOTAL',
            'harvest': total_harvest,
            'value': total_harvest * args['unit_price'] * args['frac_post_process'],
        })

        harvest_value_df = pandas.DataFrame.from_dict(harvest_value)
        # Rearrange the columns
        harvest_value_df = harvest_value_df[['subregion', 'harvest', 'value']]

        LOGGER.info('HARVEST VALUE price=%s, frac_post_process=%s',
            args['unit_price'], args['frac_post_process'])
        print(harvest_value_df)
    else:
        LOGGER.info('Valuation disabled')

    # Totals-per-timestep
    # produce a dataframe of (timestep, spawners, recruits, harvest, value, equilibrium)
    # If no valuation, leave out the valuation column.
    df_spawners_recruits = pandas.DataFrame.from_dict({
        'Timestep': list(range(0, n_timesteps + 1)),
        'Total spawners': total_spawners[:-1],
        'Total recruits': total_recruits[:-1]})
    df_spawners_recruits = df_spawners_recruits[['Timestep', 'Total spawners', 'Total recruits']]
    total_harvest_series = harvest_df.sum(axis=1)
    total_harvest_df = pandas.DataFrame({'Total harvest': total_harvest_series})
    df_spawners_recruits = df_spawners_recruits.join(total_harvest_df)

    # If we're doing valuation, add in a total harvest value column.
    if args['val_cont']:
        harvest_value_df = pandas.DataFrame({
            'Harvest value': total_harvest_series * args['unit_price'] * args['frac_post_process']
        })
        df_spawners_recruits = df_spawners_recruits.join(harvest_value_df)

    # Should this always be based on value?  Or can this be based on weight?
    equilibrated_series = [None] * 10
    is_equilibrated = [None] * 10
    for timestep in range(10, n_timesteps+1):
        equilibration = (
            total_harvest_series[timestep-9:timestep+1].mean() /
            total_harvest_series[timestep])
        equilibrated_series.append(equilibration)
        if numpy.isclose(equilibration, 1.0, rtol=1e-4):
            is_equilibrated.append('Y')
        else:
            is_equilibrated.append('N')
    equilibrated_df = pandas.DataFrame({
        'Equilibration checker': equilibrated_series,
        'Is equilibrated': is_equilibrated,
    })
    df_spawners_recruits = df_spawners_recruits.join(equilibrated_df)

    print(df_spawners_recruits)

    return tuple(df_spawners_recruits[['Total spawners', 'Total harvest']].iloc[[-1]].values[0])



if __name__ == '__main__':
    #lobster_jess()
    lobster()
    shrimp()
    dungeness_crab()
    blue_crab()
