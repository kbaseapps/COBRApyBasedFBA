import os
import json
import math
import jinja2
import numpy as np

# Helper function
def round_float_str(s, precision=6):
    """
    Given a string containing floats delimited by whitespace
    round each float in the string and return the rounded string.

    Parameters
    ----------
    s : str
        String containing floats
    precision : int
        Number of decimals to round each float to
    """
    round_str = ''
    for token in s.split():
        try:
            f = round(float(token), precision)
            round_str += str(f)
        except:
            round_str += token
        round_str += ' '
    return round_str

# Helper functions for formating
yes_no_format = lambda x: 'Yes' if x else 'No'
missing_format = lambda x: x if x else '-'
nan_format = lambda x: '-' if np.isnan(x) else x
round_format = lambda x: nan_format(round(x, 6))
equation_format = lambda rct: round_float_str(rct.build_reaction_string(use_metabolite_names=True))

# Helper function to determine reaction class
def class_formater(rct_id, fva_sol):
    min_zero = math.isclose(fva_sol.minimum[rct_id], 0, abs_tol=1e-07)
    max_zero = math.isclose(fva_sol.maximum[rct_id], 0, abs_tol=1e-07)

    min_ = 0 if min_zero else fva_sol.minimum[rct_id]
    max_ = 0 if max_zero else fva_sol.maximum[rct_id]

    if min_zero and max_zero:
        return 'blocked'
    if min_ > 0 or max_ < 0:
        return 'essential'
    #if min_ < 0 or max_ > 0:
    return 'functional'

# Helper function to format reaction/ex-reaction display data
def reaction_formater(fba_sol, fva_sol, ex):
    """ex specifies exchange reaction"""
    if fva_sol is None:
        return json.dumps([])

    # Select either exchange or normal reactions
    if ex:
        rct_ids = fva_sol.loc[fva_sol.index.str[:2] == 'EX'].index
    else:
        rct_ids = fva_sol.loc[fva_sol.index.str[:2] != 'EX'].index

    # Get reaction objects from ids
    rcts = map(model.reactions.get_by_id, rct_ids)

    return json.dumps([{'id': rct_id,
                        'flux': round_format(fba_sol.fluxes[rct_id]),
                        'min_flux': round_format(fva_sol.minimum[rct_id]),
                        'max_flux': round_format(fva_sol.maximum[rct_id]),
                        'class': class_formater(rct_id, fva_sol),
                        'equation': equation_format(rct),
                        'name': missing_format(rct.name)}
                       for rct_id, rct in zip(rct_ids, rcts)])

# Helper function for formatting model summary
def model_summary(model):

    df = model.summary().to_frame()

    def rct_name(rct_id):
        if rct_id is np.nan:
            return '-'

        try:
            name = model.reactions.get_by_id('EX_' + rct_id).name
        except:
            pass
        else:
            return name + f'\n({rct_id})'

        try:
            name = model.reactions.get_by_id(rct_id).name
        except:
            pass
        else:
            return name + f'\n({rct_id})'

        try:
            name = model.metabolites.get_by_id(rct_id).name
        except:
            pass
        else:
            return name + f'\n({rct_id})'

        return '-'

    df[('IN_FLUXES',  'ID')] = df[('IN_FLUXES',    'ID')].apply(rct_name)
    df[('IN_FLUXES',  'FLUX')] = df[('IN_FLUXES',  'FLUX')].apply(round_format)
    df[('OUT_FLUXES', 'ID')] = df[('OUT_FLUXES',   'ID')].apply(rct_name)
    df[('OUT_FLUXES', 'FLUX')] = df[('OUT_FLUXES', 'FLUX')].apply(round_format)
    df[('OBJECTIVES', 'ID')] = df[('OBJECTIVES',   'ID')].apply(rct_name)
    df[('OBJECTIVES', 'FLUX')] = df[('OBJECTIVES', 'FLUX')].apply(round_format)

    return [row[1:] for row in df.itertuples()]

# Helper function for formatting ATP summary
def atp_summary_formatter(model):
    """Returns list of ATP summary values if metabolites are found
       or a message stating they could not be found. Also return a
       bool specicifying whether or not summary exists."""
    # Select ATP metabolite
    if 'atp_c' in model.metabolites:
        df = model.metabolites.atp_c.summary().to_frame()
    elif 'ATP_c' in model.metabolites:
        df = model.metabolites.ATP_c.summary().to_frame()
    elif 'cpd00002_c0' in model.metabolites:
        df = model.metabolites.cpd00002_c0.summary().to_frame()
    else:
        # Empty ATP summary
        msg = 'Could not find atp_c, ATP_c or cpd00002_c0 in metabolites. ' \
              'Add either of these metabolites to the model in ' \
              'order to display an ATP summary.'
        return msg, False

    print(df)

    rcts = [(model.reactions.get_by_id(rct_id), rct_id)
            for rct_id in df.index.get_level_values(1)]

    df['FLUX'] = df['FLUX'].apply(round_format)
    df['PERCENT'] = df['PERCENT'].apply(lambda x: f'{round(x, 2)}%')
    df['SIDE']= df.index.get_level_values(0)
    df['NAME_ID'] = [rct.name + f'\n({rct_id})' for rct, rct_id in rcts]
    df['REACTION_STRING'] = [equation_format(rct) for rct, _ in rcts]

    atp_summary = [list(row) for row in zip(df['SIDE'], df['NAME_ID'], df['PERCENT'],
                                            df['FLUX'], df['REACTION_STRING'])]
    return atp_summary, True

# Helper function for formating essential genes
def essential_genes_formatter(model, essential_genes):
    if not essential_genes:
        return json.dumps([])

    return json.dumps([{'name': missing_format(gene.id),
                        'essential': yes_no_format(gene in essential_genes)}
                      for gene in model.genes])

# Call this function to build the report
def build_report(pipeline, model, fba_sol, fva_sol,
                 essential_genes, model_id, media_id):
    """Build output report and return string of html."""

    atp_summary, is_atp_summary = atp_summary_formatter(model)

    # Formating for objective value
    obj_value = round_format(fba_sol.fluxes[pipeline.target_reaction])
    obj_units = 'gm/gm CDW hr' if 'biomass' in pipeline.target_reaction else 'mmol/gm CDW hr'

    context = {'summary':     model_summary(model),
               'atp_summary': {'is_atp_summary': is_atp_summary, 'summary': atp_summary},
               'overview':    [{'name': 'Model',                    'value': model_id},
                               {'name': 'Media',                    'value': media_id},
                               {'name': 'Optimization status',      'value': fba_sol.status},
                               {'name': 'Objective',                'value': model.objective},
                               {'name': 'Target objective value',   'value': f'{obj_value} ({obj_units})'},
                               {'name': 'Number of reactions',      'value': len(model.reactions)},
                               {'name': 'Number of compounds',      'value': len(model.metabolites)},
                               {'name': 'FBA type',                 'value': pipeline.fba_type},
                               {'name': 'FVA type',                 'value': pipeline.fva_type},
                               {'name': 'FVA fraction of optimum',  'value': pipeline.fraction_of_optimum_fva},
                               {'name': 'All reversible reactions', 'value': yes_no_format(pipeline.is_all_reversible)},
                               {'name': 'Single gene KO',           'value': yes_no_format(pipeline.is_single_ko)},
                               {'name': 'Gene KO',                  'value': len(pipeline.feature_ko_list)},
                               {'name': 'Reaction KO',              'value': len(pipeline.reaction_ko_list)},
                               {'name': 'Custom bounds',            'value': len(pipeline.custom_bound_list)},
                               {'name': 'Media supplement',         'value': len(pipeline.media_supplement_list)},
                               {'name': 'Solver',                   'value': pipeline.solver.upper()}
                               ],
               'reaction_tab': {
                   'is_reactions': fva_sol is not None,
                   'reactions': reaction_formater(fba_sol, fva_sol, ex=False),
                   'help': 'Select FVA setting and rerun to produce results.'
                },
               'ex_reaction_tab': {
                   'is_reactions': fva_sol is not None,
                   'reactions': reaction_formater(fba_sol, fva_sol, ex=True),
                   'help': 'Select FVA setting and rerun to produce results.'
                },
                'essential_genes_tab': {
                    'is_essential_genes': len(essential_genes) > 0,
                    'essential_genes': essential_genes_formatter(model, essential_genes),
                    'help': 'Select simulate all single KO to produce results.'
                }
           }

    template_file = 'template.html'
    # Directory this file is in
    template_dir = os.path.dirname(os.path.realpath(__file__))

    env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(template_dir),
            autoescape=jinja2.select_autoescape(['html', 'xml']))

    # Return string of html
    return env.get_template(template_file).render(context)
