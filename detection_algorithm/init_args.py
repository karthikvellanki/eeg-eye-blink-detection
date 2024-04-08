def args_init(delta_uV):
    args = {}
    args['mintab'], args['maxtab'] = [], []
    args['mn'], args['mx'] = float("inf"), -1*float("inf")
    args['mnpos'], args['mxpos'] = None, None
    args['min_left'], args['min_right'] = [], []
    args['lookformax'] = True
    args['delta'] = delta_uV
    return args