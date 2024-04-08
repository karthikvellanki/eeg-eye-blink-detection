def identify_signal_peaks(time, value, args):
    foundMin = False
    if value > args['mx']:
        args['mx'] = value
        args['mxpos'] = time
    if value < args['mn']:
        args['mn'] = value
        args['mnpos'] = time
    if args['lookformax']:
        if value < args['mx'] - args['delta']:
            args['maxtab'].append([args['mxpos'], args['mx']])
            args['mn'] = value
            args['mnpos'] = time
            args['lookformax'] = False
    else:
        if value > args['mn'] + args['delta']:
            args['mintab'].append([args['mnpos'], args['mn']])
            args['min_left'].append([-1, -1])
            args['min_right'].append([-1, -1])
            args['mx'] = value
            args['mxpos'] = time
            args['lookformax'] = True                
            foundMin = True
    return foundMin