#!/usr/bin/env python
#
# get_all_obj.py
"""

diagnostic tool to search for memory leaks.

Gathers all objects not found by the Garbage Collector
and counts them by type.

"""
import gc

# Recursively expand slist's objects
# into olist, using seen to track
# already processed objects.
def _getr(slist, olist, seen):
  for e in slist:
    if id(e) in seen:
      continue
    seen[id(e)] = None
    olist.append(e)
    tl = gc.get_referents(e)
    if tl:
      _getr(tl, olist, seen)

# The public function.
def get_all_objects():
  """Return a list of all live Python
  objects, not including the list itself."""
  gc.collect()
  gcl = gc.get_objects()
  olist = []
  seen = {}
  # Just in case:
  seen[id(gcl)] = None
  seen[id(olist)] = None
  seen[id(seen)] = None
  # _getr does the real work.
  _getr(gcl, olist, seen)
  return olist

def print_all_objects():
    tp  = {}
    obj = get_all_objects()
    for o in obj:
        stype = str(type(o))
        if stype == "<type 'instance'>":
            if hasattr(o, '__module__'):
                stype = 'from module:'+str(o.__module__)
            else:
                stype = str(o)
        tp.setdefault(stype, 0)
        tp[stype] += 1
    rr = zip(tp.values(), tp.keys())
    rr.sort()
    for a, b in rr:
        print a, b
    print 'total', len(obj)
        
