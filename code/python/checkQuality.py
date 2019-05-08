
import pandas
import time


# Identify gapfilled reactions
def findGapFilled(model, exclude=[]):
    gapfilled = []
    if not type(exclude) is list:
        exclude = [exclude]
        
    for index in model.reactions:
        if len(list(index.genes)) == 0:
            if not index in model.boundary:
                if not index.id in exclude:
                    gapfilled.append(index.id)
    
    if len(gapfilled) > 0:
        print(str(len(gapfilled)) + ' reactions not associated with genes')
    
    return gapfilled


# Checks which metabolites are generated for free
def checkFreeMass(raw_model):
    
    with raw_model as model:
        
        for index in model.boundary:
            model.reactions.get_by_id(index.id).lower_bound = 0.
              
        demand_metabolites = [x.reactants[0].id for x in model.demands if len(x.reactants) > 0] + [x.products[0].id for x in model.demands if len(x.products) > 0]

        free = []
        for index in model.metabolites: 
            if index.id in demand_metabolites:
                continue
            else:
                demand = model.add_boundary(index, type='demand')
                model.objective = demand
                obj_val = model.slim_optimize(error_value=0.)
                if obj_val > 1e-8:
                    free.append([index.id, obj_val])
                model.remove_reactions([demand])
    
    if len(free) > 0:
        print(str(len(free)) + ' metabolites are generated for free')
 
    return(free)


# Check for mass and charge balance in reactions
def checkBalance(raw_model, exclude=[]):
    
    with raw_model as model:
        imbalanced = []
        mass_imbal = 0
        charge_imbal = 0
        elem_set = set()
        for metabolite in model.metabolites:
            try:
                elem_set |= set(metabolite.elements.keys())
            except:
                pass
        
        if not type(exclude) is list: exclude = [exclude]
        for index in model.reactions:
            if index in model.boundary or index.id in exclude:
                continue
                     
            else:
                try:
                    test = index.check_mass_balance()
                except ValueError:
                    continue
                    
                if len(list(test)) > 0:
                    imbalanced.append(index.id)
                    
                    if 'charge' in test.keys():
                        charge_imbal += 1
                    if len(set(test.keys()).intersection(elem_set)) > 0:
                        mass_imbal += 1
                
    if mass_imbal != 0:
        print(str(mass_imbal) + ' reactions are mass imbalanced')
    if charge_imbal != 0:
        print(str(charge_imbal) + ' reactions are charge imbalanced')
    
    return(imbalanced)


# Identifies blocked reactions, 1% cutoff for fraction of optimum
def blockedReactions(model):
    
    fva = flux_variability_analysis(model, fraction_of_optimum=0.01)
    noflux = (fva["maximum"].abs() < 1e-8) & (fva["minimum"].abs() < 1e-8)    
    blocked = noflux[noflux==True].index.tolist()

    if noflux.sum() != 0:
        print(str(noflux.sum()) + ' reactions are blocked')
        
    return blocked


# Checks the quality of models by a couple metrics and returns problems
def checkQuality(model, exclude=[]):
    
    start_time = time.time()
    
    if model.name != None:
        model_name = model.name
    else:
        model_name = 'model'
    
    gaps = findGapFilled(model, exclude)
    freemass = checkFreeMass(model)
    balance = checkBalance(model, exclude)
    blocked = blockedReactions(model)
    
    test = gaps + freemass + balance + blocked
    if len(test) == 0:
        print('No inconsistencies detected')
    
    duration = int(round(time.time() - start_time))
    print('Took ' + str(duration) + ' seconds to analyze ' + model_name) 
    
    return gaps, freemass, balance, blocked
