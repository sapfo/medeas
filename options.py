import options

options.TESTING = True  # Shows some debugging info and many plots
options.FST = False  # Also calculates F_ST
options.BOOTRUNS = 10  # How many bootstrap runs we need.
options.SIMULATION = False  # Are we running on simulated datas
options.K_OVERRIDE = 6 #number of population shoule be overwritten? k == 0 -> no; K > 0 -> number of population
options.MATHPLOTLIB_BACKEND = 'agg' #which backend should be used by mathplotlib, Might be "agg" or ???macosx???
options.VERBOSE = 2 #0 only minimal output, 1 important step, 2 much more stuff
