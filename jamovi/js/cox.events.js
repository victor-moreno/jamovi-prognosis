// Interaction builder for the Cox analysis.
// The pool of variables that can be crossed is the union of the Factors and
// Covariates boxes; terms that reference a removed variable are pruned.
// Adapted from jsurvival (GPL >= 2) multisurvival.events.js, jus 2.0 API.

const events = {
    update: function(ui) {
        calcInteractionPool(ui, this);
    },

    onChange_predictors: function(ui) {
        calcInteractionPool(ui, this);
    },

    onUpdate_interactionSupplier: function(ui) {
        let vars = collectPredictors(ui, this);
        ui.interactionSupplier.setValue(this.valuesToItems(vars, FormatDef.variable));
    }
};

let collectPredictors = function(ui, context) {
    let a = context.cloneArray(ui.factors.value(), []);
    let b = context.cloneArray(ui.covs.value(), []);
    return a.concat(b);
};

let calcInteractionPool = function(ui, context) {
    let vars = collectPredictors(ui, context);
    ui.interactionSupplier.setValue(context.valuesToItems(vars, FormatDef.variable));

    let varsDiff = context.findChanges("predictorList", vars, true, FormatDef.variable);
    let termsList = context.cloneArray(ui.interactions.value(), []);
    let changed = false;
    for (let i = 0; i < varsDiff.removed.length; i++) {
        for (let j = 0; j < termsList.length; j++) {
            if (FormatDef.term.contains(termsList[j], varsDiff.removed[i])) {
                termsList.splice(j, 1);
                changed = true;
                j -= 1;
            }
        }
    }
    if (changed)
        ui.interactions.setValue(termsList);
};

module.exports = events;
