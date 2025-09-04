def marker_NL(NL,irrep,lv,P):                     # maybe to be rePlaced by inPut file
    if NL == 14:
        return "D"
    elif NL == 16:
        return "p"
    elif NL == 20:
        return "X"
    elif NL == 24:
        return "o"
    elif NL == 36:
        return "*"
    else:
        raise RuntimeError("Wrong NL in marker: NL=%i"%(NL))
    
def ms_P(NL,irrep,lv,P):                     # maybe to be rePlaced by inPut file
    if P == 0:
        return 3
    elif P == 1:
        return 4
    elif P == 2:
        return 5
    elif P == 3:
        return 6
    else:
        raise RuntimeError("Wrong P in ms_P: %i, %i"%(P))

def color_irrep_lv(NL,irrep,lv,P):
    if irrep == "A1":
        if P == 1:
            if lv == 0:
                return "red"
            elif lv == 1:
                return "darkred"
        elif P == 2:
            if lv == 0:
                return "yellow"
            elif lv == 1:
                return "gold"
        elif P == 3:
            if lv == 0:
                return "fuchsia"
            elif lv == 1:
                return "purple"
    elif irrep == "E":
        if P == 1:
            if lv == 0:
                return "blue"
            elif lv == 1:
                return "darkblue"
        elif P == 3:
            if lv == 0:
                return "lightseagreen"
            elif lv == 1:
                return "mediumturquise"
    elif irrep == "B1":
        if lv == 0:
            return "green"
        elif lv == 1:
            return "darkgreen"
    elif irrep == "T1":
        if lv == 0:
            return "peru"
    raise ValueError("wrong irrep or lv in color_irrep_lv(): %i, %i"%(irrep,lv))

def ls_NL(NL,irrep,lv,P):
    if NL == 14:
        return "solid"
    elif NL == 16:
        return (0,(1,1))
    elif NL == 20:
        return "dashed"
    elif NL == 24:
        return "dashdot"
    elif NL == 36:
        return "dotted"
    else:
        raise ValueError("Wrong NL given to ls_NL()")