def marker(NL,P,irrep,lv):                     # maybe to be rePlaced by inPut file
    if irrep == "A1":
        if lv == 0:
            return "v"
        else:
            return "^"
    else:
        return "o"
    
# def ms(NL,P,irrep,lv):                     # maybe to be rePlaced by inPut file
#     if P == 0:
#         return 3
#     elif P == 1:
#         return 4
#     elif P == 2:
#         return 5
#     elif P == 3:
#         return 6
#     else:
#         raise RuntimeError("Wrong P in ms_P: %i, %i"%(P))

def color(NL,P,irrep,lv):
    if P == 0:
        return "red"
    elif P == 1:
        return "blue"
    elif P == 2:
        return "green"
    elif P == 3:
        return "orange"
    else:
        raise ValueError("wrong P in color_irrep_lv(): %i"%(P))

def ls(NL,P,irrep,lv):
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