

import uedge
from uedge import uedge_lists as ul
import hashlib

def findvar(var):
    for p in ul.packages:
        if var in p.varlist():
            return p,p.getpyobject(var)
    return None,None

class watch:
        def __init__(self,var):
           pkg,pvar = findvar(var)
           if pvar is not None:
              try:
                 m=hashlib.md5()
                 self.sval = str(pvar)
                 m.update(self.sval.encode('UTF-8'))
                 self.hash = m.hexdigest()
                 self.val = pvar.copy()
              except:
                 self.val = pvar
              self.var = var
           else:
              #print("No variable in Uedge found called ",var)
              self.var = None
        def changed(self):
           pkg,pvar = findvar(self.var)
           if pvar is not None:
              try:
                  m=hashlib.md5()
                  self.sval = str(pvar)
                  m.update(self.sval.encode('UTF-8'))
                  hash = m.hexdigest()
                  if hash != self.hash: return True
                  if (self.val != pvar).any(): return True
              except:
                  if self.val != pvar: return True
           return False 


class watchall:
     def __init__(self):
         self.watching = []
         for p in ul.packages:
             for v in p.varlist():
                 self.watching.append(watch(v))
     def printchanged(self):
         for w in self.watching:
              if w.changed():
                  print(w.var)
     def check(self):
         changed = []
         for w in self.watching:
              if w.changed():
                  changed.append(w.var)
         return changed



       
