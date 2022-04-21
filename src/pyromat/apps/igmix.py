

def igmix(sid, contents, bymass=True, doc=''):
    """Create an ideal gas mixture
    
    igm = igmix(sid, contents)
        OR
    igm = igmix(sid, contents, bymass=False)
    
Creates an igmix instance with a species id specified by sid and contents
specified by contents.  The result will be returned, but it will not be
added to the PYroMat data system, so it will not be visible by get() or
info(), and it will not be available if PYroMat is reloaded.

To make the new mixture a permanent addition to the PYroMat data system,
see the new() and updatefiles() functions in the dat module.

** Species ID **
If the sid does not begin with "ig." to signifiy that it belongs to the
ideal gas collection, it will be automatically added.  If the substance
is to be added to the PYroMat data system, its SID must be unique.  
Otherwise, the sid will be ignored.


"""
    # Force this into the ideal gas collection
    if not sid.startswith('ig.'):
        sid = 'ig.' + sid

    # Verify that igmix is defined
    if 'igmix' not in pm.reg.registry:
        raise pm.utility.PMDataError('The igmix class does not seem to be available in the current registry. Aborting.')

    # Verify the contents are all valid ideal gas substances
    for ss in contents:
        if ss not in pm.dat.data or not pm.dat.data[ss].data['id'].startswith('ig.'):
            message = 'The igmix constituent was not a recognized ideal gas: ' + ss
            pm.utility.print_error(message)
            raise pm.utility.PMDataError(message)

    # Build the data dictionary to simulate JSON data from a file
    data = {'id':sid, 'contents':contents, 'bymass':bymass, 'class':'igmix', 'doc':doc, 'fromfile':None}

    # Generate the igmix instance and return
    return pm.reg.registry['igmix'](data)




class digmix:
    """Placeholder for the future Dynamic Ideal Gas Mixture class
    
This is a planned class for calculating the properties of a variable 
mixture of ideal gases."""
    pass
