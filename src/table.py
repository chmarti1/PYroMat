#!../../bin/python
import pyromat as pm

s = list(pm.search(collection='ig'))
s.sort(key=lambda x:x.sid())

with open('table', 'w') as fd:
    for this in s:
        names = this.names()
        if names:
            fd.write(f'<tr><td><span class="code">{this.sid()}</span></td><td>{names[0]}</td></tr>\n')
        else:
            sid = this.sid()
            name = sid.split('.')[1]
            fd.write(f'<tr><td><span class="code">{sid}</span></td><td>{name}</td></tr>\n')

