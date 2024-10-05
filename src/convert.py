#!/usr/bin/python3

import json
import os,sys

def convert1(original):
    output = []
    # Loop over each polynomial coefficient group
    for pgroup in original:
        prex = pgroup[0]
        postx = pgroup[1]
        pcoef = {'coef':pgroup[2:]}
        if prex!=1:
            pcoef['pre'] = prex
        if postx!=0:
            pcoef['post'] = postx
        output.append(pcoef)
        
    if len(output) == 1:
        return output[0]
    return output

def convert2(original):
    output = []
    # Loop over each polynomial coefficient group
    for pgroup in original:
        prex,prey = pgroup[0]
        postx,posty = pgroup[1]
        pcoef = {'coef':pgroup[2:]}
        if prex!=1 or prey!=1:
            pcoef['pre'] = [prex, prey]
        if postx!=0 or posty!=0:
            pcoef['post'] = [postx, posty]
        output.append(pcoef)
        
    if len(output) == 1:
        return output[0]
    return output

def move(src, dest, name, dname=None, conv=None):
    """Move a data element from one dictionary to another
    move(src, dest, name,...)

src     Source dictionary
dest    Destination dictionary
name    Name of the element in the dictionary
--- Optional Params ---
dname   Name of the element in the destination dictionary (change name)
conv    Conversion function to modify the data prior to moving
--- ---

If the name is not found in the source dictionary, the move will not 
happen.
"""
    if dname is None:
        dname = name
    data = src.get(name)
    if data is not None:
        if conv is not None:
            dest[dname] = conv(data)
        else:
            dest[dname] = data

if __name__ == '__main__':
    filename = sys.argv[1]
    
    pat,fn = os.path.split(filename)
    newfilename = os.path.join(pat, 'new_' + fn)
    
    with open(filename,'r') as fd:
        data = json.load(fd)

    # Modify data... ALL polynomials are affected
    
    # AOgroup 
    #   Convert coef0 to group0
    #   Rename coef1 to group1
    group = data.pop('AOgroup')
    data['IGgroup'] = {
        'Tscale':group['Tscale'],
        'dscale':group['dscale'],
    }
    move(group, data['IGgroup'], 'logt')
    move(group, data['IGgroup'], 'tlogt')
    move(group, data['IGgroup'], 'coef0', 'group0', convert1)
    move(group, data['IGgroup'], 'coef1', 'group1')
    
    # ARgroup
    #   Convert each element of coef0 to group0
    group = data.pop('ARgroup')
    data['Rgroup'] = {
        'Tscale':group['Tscale'],
        'dscale':group['dscale'],
        'group0':[],
    }
    for this in group['coef0']:
        data['Rgroup']['group0'].append(convert2(this))
    move(group, data['Rgroup'], 'coef1', 'group1')
    move(group, data['Rgroup'], 'coef2', 'group2')
    
    # DSLgroup
    #   coef
    data['DSLgroup']['poly'] = convert1(data['DSLgroup'].pop('coef'))
    
    # DSVgroup
    #   coef
    data['DSVgroup']['poly'] = convert1(data['DSVgroup'].pop('coef'))
    
    # PSgroup
    #   coef
    data['PSgroup']['poly'] = convert1(data['PSgroup'].pop('coef'))

    with open(newfilename,'w') as fd:
        json.dump(data, fd, indent=4)
