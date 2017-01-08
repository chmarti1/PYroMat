import pyromat as pyro
import re

steam = pyro.get('steam')


def add_tag(tag,iscale=1,jscale=1):
    datafile = tag + '.dat'

    fixer = re.compile('\\xe2\\x80\\x93')

    C = []
    with open(datafile,'r') as ff:
        for thisline in ff:
            thisline = fixer.sub('-',thisline)
            temp = thisline.split()
            while len(temp)>3:
                temp[2]+=temp.pop(3)
            C.append([int(iscale*float(temp[0])), int(jscale*float(temp[1])), float(temp[2])])

    steam.data[tag] = C



add_tag('th1')
add_tag('ts1')
add_tag('th2a')
add_tag('th2b')
add_tag('th2c')
add_tag('ts2a',iscale=4)
add_tag('ts2b')
add_tag('ts2c')
steam.data['b2bc'] = [0.90584278514723e3, -0.67955786399241,
    0.12809002730136e-3, 0.26526571908428e4, 0.45257578905948e1]


