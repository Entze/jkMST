const loglevels = [Logging.Debug, Logging.Info, Logging.Warn, Logging.Error]

@enum Solver glpk scip cplex

@enum SolvingMode mtz scf mcf cec dcc
