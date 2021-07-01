#!/usr/bin/python3

# Add atoms to all species
import pyromat as pm
import os,sys

with open('addcontents.log', 'w') as ff:

    for spid,species in pm.dat.data.items():
        ff.write(spid + '    ')
        # Ignore IGMIX classes
        if isinstance(species, pm.reg.registry['igmix']):
            ff.write('Ignoring\n')
        elif 'atoms' not in species.data:
            collection,name = spid.split('.')
            # Construct the chemical formula from the name portion of the id string
            atoms = {}
            atom = None
            qty = ''
            sign = 1
            # Create a special case for the free electron
            if name.startswith('e'):
                atoms = {'e':1}
            else:
                for char in name:
                    # If this is the beginning of a new atom
                    if 'A' <= char <= 'Z':
                        # If this is not the first atom, stow the results of the last
                        if atom is not None:
                            if qty:
                                atoms[atom] = int(qty) * sign
                            else:
                                atoms[atom] = sign
                        atom = char
                        qty = ''
                        sign = 1
                    # If the species is an ion
                    elif char == '+' or char == '-':
                        # If this is not the first atom, stow the results of the last
                        if atom is not None:
                            if qty:
                                atoms[atom] = int(qty) * sign
                            else:
                                atoms[atom] = 1 * sign
                        atom = 'e'
                        qty = ''
                        sign = 1 if char == '-' else -1

                    # If this is a two-letter atom
                    elif 'a' <= char <= 'z':
                        atom = atom + char
                    elif '0' <= char <= '9':
                        qty = qty + char
                # Finally, stow the last atom
                if qty:
                    atoms[atom] = int(qty) * sign
                else:
                    atoms[atom] = sign
                
                for atom in atoms:
                    ff.write(atom + ':' + str(atoms[atom]) + '  ')
                ff.write('\n')
                
            species.data['atoms'] = atoms
        else:
            ff.write('\n')
            
pm.dat.updatefiles()
