
# Latitude plotting fuction to get N/S labels right:
# see https://jakevdp.github.io/PythonDataScienceHandbook/04.10-customizing-ticks.html
#

def NS_format_func(value, tick_number):
    
    # add north/south latitude labels

    N = int(value)
    if N < 0:
        if N == -90:
            return ""
        else:
            return r"${0}S$".format(-N)
    elif N == 0:
        return "0"
    elif N > 0 :
        return r"${0}N$".format(N)
    else:
        return "oops!"

# Pressure axis major label formatter function to get labels right:

def millibar_format_func(value, tick_number):
    N = int(value)
    if N == 1000:
        return ""
    else:
        return r"${0}$".format(N)
