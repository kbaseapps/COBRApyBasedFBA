import math
import cobra
import cobrakbase
import numpy as np
import pandas as pd
from cobrakbase.modelseed.utils import atom_count
from cobrakbase.core.kbase_fba_builder import KBaseFBABuilder

class FBAPipeline:

    # Bound for making reactions reversible
    MAX_BOUND = 1000

    # Tuple of uptake atoms
    UPTAKE_ATOMS = ('C', 'N', 'P', 'S', 'O')
    
    def __init__(self):
        # If true, run FVA algorithm on model.
        self.is_run_fva = False

        # If true, make all reactions in model reversible.
        self.is_all_reversible = False

        # If true, run gene knockout algorithm on each gene
        # to determine the set of the essential genes.
        self.is_single_ko = False

        # Sets FBA type, either 'pFBA', 'loopless FBA' or 'FBA'
        self.fba_type = 'pFBA'

        # If true, run loopless FVA with CycleFreeFlux.
        self.is_loopless_fva = False

        # If true, the objective will be minimized.
        # Otherwise the objective is maximized.
        self.is_minimize_objective = False

        # For pFBA, Must be <= 1.0. Requires that the objective value
        # is at least the fraction times maximum objective value.
        self.fraction_of_optimum_pfba = 1.0

        # For FVA, Must be <= 1.0. Requires that the objective value
        # is at least the fraction times maximum objective value.
        self.fraction_of_optimum_fva = 0.1

        # Specify optimization solver, options: [coinor_cbc, glpk].
        self.solver = 'coinor_cbc'

        # Specify objective reaction to optimize.
        self.target_reaction = ''

        # Max uptakes for exchange reactions
        self.max_uptakes = {atom: None for atom in self.UPTAKE_ATOMS}

        # Used with complete media. Default is 0, if complete media
        # change from 0 to 100. Otherwise if user specifies
        # default_max_uptake != 0, every exchange flux has a lower
        # bound of negative default_max_uptake.
        self.default_max_uptake = 0.

        # Knocks out specified genes (which could knockout reactions)
        self.feature_ko_list = []

        # Knockout specified reactions.
        self.reaction_ko_list = []

        # Custom bounds to add to reaction in model
        self.custom_bound_list = []

        # TODO: write comment. should this be removed?
        self.media_supplement_list = []

        # Kbase ID of the returning FBA object.
        self.output_id = ''

        # Kbase workspace name where FBA app is run.
        self.workspace = ''

        # NOTE: boolean variables are passed in as ints
        #       with 0,1 value during runtime.

    @staticmethod
    def fromKBaseParams(params):
        p = FBAPipeline()

        # Kbase environment params
        p.output_id = params['fba_output_id']
        p.workspace = params['fbamodel_workspace']

        # Optimization and algorithmic params
        p.solver = params['solver']
        p.is_run_fva = params['fva']
        p.fba_type = params['fba_type']
        p.is_single_ko = params['simulate_ko']
        p.is_loopless_fva = params['loopless_fva']
        p.target_reaction = params['target_reaction']
        p.is_all_reversible = params['all_reversible']
        p.is_minimize_objective = params['minimize_objective']
        p.fraction_of_optimum_fva = params['fraction_of_optimum_fva']
        p.fraction_of_optimum_pfba = params['fraction_of_optimum_pfba']

        # Uptakes
        p.max_uptakes['C'] = params['max_c_uptake']
        p.max_uptakes['N'] = params['max_n_uptake']
        p.max_uptakes['P'] = params['max_p_uptake']
        p.max_uptakes['S'] = params['max_s_uptake']
        p.max_uptakes['O'] = params['max_o_uptake']
        p.default_max_uptake = params['default_max_uptake']

        # Check if list params contain data. If so parse, else use default []
        if params['media_supplement_list']:
            p.media_supplement_list = params['media_supplement_list'].split(',')
        if params['reaction_ko_list']:
            p.reaction_ko_list = params['reaction_ko_list'].split(',')
        if params['feature_ko_list']:
            p.feature_ko_list = params['feature_ko_list'].split(',')
        # custom_bound_list input format [{'custom_reaction_id': ['rct1'], 'custom_lb': float, 'custom_ub': float}]
        if params['custom_bound_list']:
            p.custom_bound_list = [(item['custom_reaction_id'][0], item['custom_lb'], item['custom_ub'])
                                    for item in params['custom_bound_list']]

        return p

    def configure_media(self, model, media):

        if media.name == 'Complete':
            self.default_max_uptake = 100.

        if not math.isclose(self.default_max_uptake, 0):
            for ex_flux in model.medium:
                model.reactions.get_by_id(ex_flux).lower_bound = -1 * self.default_max_uptake

        # Only add constraints when user specifies to, otherwise max_uptakes will be None
        constrs = {atom: 0. for atom in self.UPTAKE_ATOMS if self.max_uptakes[atom] is not None}

        # Return early if user doesn't specifies any max uptakes
        if not constrs:
            return

        for rct_id in model.medium:

            ex_rct = model.reactions.get_by_id(rct_id)
            compound = list(ex_rct.metabolites)[0]
            cmp_atoms = atom_count(compound.formula)

            for atom in constrs:
                atom_occurences = cmp_atoms.get(atom)
                if atom_occurences:
                    constrs[atom] += ex_rct.reverse_variable * atom_occurences

        for atom in constrs:
            model.add_cons_vars(
                model.problem.Constraint(constrs[atom],
                                         lb=0,
                                         ub=self.max_uptakes[atom]))

    def run(self, model, media):
        """This function mutates model."""

        # Select optimization solver
        model.solver = self.solver
        if self.solver == 'coinor_cbc':
            # Use all available processors
            model.solver.configuration.threads = -1

        # If specified, make all reactions reversible
        if self.is_all_reversible:
            for rct in model.reactions:
                rct.upper_bound = self.MAX_BOUND
                rct.lower_bound = self.MAX_BOUND * -1

        # Add custom bounds to reactions
        for rct_id, lb, ub in self.custom_bound_list:
            if rct_id in model.reactions:
                rct = model.reactions.get_by_id(rct_id)
                rct.lower_bound, rct.upper_bound = lb, ub

        # Filter out user specified genes to ko that are not in model.
        self.feature_ko_list = list(filter(lambda gene: gene in model.genes,
                                           self.feature_ko_list))
        # Knockout specified genes
        if self.feature_ko_list:
            cobra.manipulation.delete_model_genes(model, self.feature_ko_list)

        for rct_id in self.reaction_ko_list:
            if rct_id in model.reactions:
                rct = model.reactions.get_by_id(rct_id)
                rct.lower_bound, rct.upper_bound = 0, 0

        # Update exchange reaction variable bounds based on user
        # specified max uptakes. Add max uptake contraints.
        self.configure_media(model, media)

        # Set objective
        if self.target_reaction:
            model.objective = self.target_reaction

        # Set direction of the objective function, maximize by default.
        if self.is_minimize_objective:
            model.objective.direction = 'min'

        # Compute FBA solution
        if self.fba_type == 'pFBA':
            from cobra.flux_analysis import pfba
            fba_sol = pfba(model,
                           fraction_of_optimum=self.fraction_of_optimum_pfba)
        elif self.fba_type == 'Loopless FBA':
            # Run CycleFreeFlux algorithm
            fba_sol = cobra.flux_analysis.loopless_solution(model)
        else:
            # Run vanilla FBA
            fba_sol = model.optimize()

        # If specified, compute FVA solution
        fva_sol = None
        if self.is_run_fva:
            from cobra.flux_analysis import flux_variability_analysis as fva
            fva_sol = fva(model,
                          loopless=self.is_loopless_fva,
                          fraction_of_optimum=self.fraction_of_optimum_fva)
        
        # If specified, simulate all single gene knockouts
        essential_genes = set()
        if self.is_single_ko:
            essential_genes = cobra.flux_analysis.variability. \
                              find_essential_genes(model, threshold=1e-11)

        # Convert COBRApy model to kbase format
        fba_builder = KBaseFBABuilder.fromCobra(self.output_id,
                                                model,
                                                fba_sol,
                                                media,
                                                self.workspace)
        
        # TODO: add fva_sol and essential_genes to this object in cobrakbase
        return fba_builder.build(), fva_sol, fba_sol, essential_genes

def build_report(pipeline, model, fba_sol, fva_sol,
                 essential_genes, model_id, media_id):
    """Build output report and return string of html."""
    import os
    import json
    import jinja2

    # Helper functions for formating
    yes_no_format = lambda x: 'Yes' if x else 'No'
    missing_format = lambda x: x if x else '-'
    nan_format = lambda x: '-' if x is np.nan else x
    round_format = lambda x: nan_format(round(x, 6))

    # Helper function to determine reaction class
    def class_formater(rct_id):
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
            josn.dumps([])

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
                            'class': class_formater(rct_id),
                            'equation': round_float_str(rct.build_reaction_string(use_metabolite_names=True)),
                            'name': missing_format(rct.name)}
                           for rct_id, rct in zip(rct_ids, rcts)])

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

        summary = []
        for row in df.itertuples():
            summary.append(row[1:])
        return summary

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

        rcts = [(model.reactions.get_by_id(rct_id), rct_id)
                for rct_id in df.index.get_level_values(1)]

        df['FLUX'] = df['FLUX'].apply(round_format)
        df['PERCENT'] = df['PERCENT'].apply(lambda x: f'{round(x, 2)}%')
        df['SIDE']= df.index.get_level_values(0)
        df['NAME_ID'] = [rct.name + f'\n({rct_id})' for rct, rct_id in rcts]
        df['REACTION_STRING'] = [round_float_str(rct.build_reaction_string(use_metabolite_names=True))
                                 for rct, _ in rcts]

        atp_summary = []
        for row in zip(df['SIDE'], df['NAME_ID'], df['PERCENT'],
                       df['FLUX'], df['REACTION_STRING']):
            atp_summary.append(list(row))

        return atp_summary, True

    def essential_genes_formatter(model, essential_genes):
        if not essential_genes:
            return json.dumps([])

        return json.dumps([{'name': missing_format(gene.id),
                            'essential': yes_no_format(gene in essential_genes)}
                          for gene in model.genes])

    atp_summary, is_atp_summary = atp_summary_formatter(model)

    # Formating for objective value
    obj_value = round_format(fba_sol.fluxes[pipeline.target_reaction])
    obj_units = 'gm/gm CDW hr' if 'biomass' in pipeline.target_reaction else 'mmol/gm CDW hr'

    # TODO: Get model_id, media_id from cobrokbase (currently they are kbase object ref)

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
                               {'name': 'Loopless FVA',             'value': yes_no_format(pipeline.is_loopless_fva)},
                               {'name': 'FVA fraction of opt.',     'value': pipeline.fraction_of_optimum_fva},
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
    template_dir = os.path.dirname(os.path.realpath(__file__))

    env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(template_dir),
            autoescape=jinja2.select_autoescape(['html', 'xml']))

    # Return string of html
    return env.get_template(template_file).render(context)
