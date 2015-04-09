# Structured Xgrow libary for reading, writing and converting Structured Xgrow (STXG) files.
# If you're just manipulating stxg files, you probably don't need this (just use yaml), unless
# you want to use the dump function, which produces somewhat cosmetically better output. The 
# library is primarily useful for converting to and from xgrow tile files.
#
# If you just want to run xgrow with stxg files, consider using xgrow-wrap.

version = '0.0.1'

import yaml
import re
from StringIO import StringIO
import datetime
import copy
#import warning

# Option names and their argument types.
keyvaloptions = [ ('block', int), ('size', int), ('rand', str), ('k', float), ('Gmc', float), ('Gse', float),
                  ('Gas', float),('Gam', float),('Gae', float),('Gfc', float),('T', float),('blast_rate_alpha', float),
                  ('blast_rate_beta', float),('blast_rate_gamma', float),('seed', str), ('update_rate', int),
                  ('tracefile', str), ('untiltiles', str), ('tmax', float),('emax', int),('smax', int),
                  ('smin', int), ('untiltilescount', str),('clean_cycles', int),('error_radius', float),
                  ('datafile', str),('arrayfile', str),('exportfile', str),('importfile', str),('min_strength', float) ]
keyopts = [x[0] for x in keyvaloptions]

# Options that need special handling.
specialoptions = ['stoic', 'nowindow', 'doubletile', 'vdoubletile', 'fission']

# Options that take true/fale values.
truefalseoptions = ['pause','wander','periodic']

# Classes for controlling flow/block style of various yaml elements. These are used in dump()
class blockseq( dict ): pass
def blockseq_rep(dumper, data):
    return dumper.represent_mapping( u'tag:yaml.org,2002:seq', data, flow_style=False )
class flowmap( dict ): pass
def flowmap_rep(dumper, data):
    return dumper.represent_mapping( u'tag:yaml.org,2002:map', data, flow_style=True )
yaml.add_representer(blockseq, blockseq_rep)
yaml.add_representer(flowmap, flowmap_rep)


def dump( tileset, *xargs, **pargs ):
    """
    stxg.dump( tileset, ... ) is a wrapper for yaml.dump that slightly prettifies
    the output taking into account what tilesets usually look like. It passes
    all arguments onto yaml.dump.
    """
    
    # Make tile and bond sections flow-style
    tileset['tiles'] = [ flowmap(x) for x in tileset['tiles'] ]
    tileset['bonds'] = [ flowmap(x) for x in tileset['bonds'] ]
    
    # If xgrowargs is there, make it block-style
    if 'xgrowargs' in tileset.keys(): 
        tileset['xgrowargs'] = blockseq(tileset['xgrowargs'])
    
    return yaml.dump(tileset, *xargs, **pargs)


def load( stream, *xargs, **pargs ):
    """
    stxg.load( stream, ... ) is a wrapper for yaml.load that currently does
    nothing else but call yaml.load
    """
    return yaml.load( stream, *xargs, **pargs )


def to_xgrow( stxg, stream=None ):
    """
    Given an stxg structure/dict, create an xgrow tile file, and return it
    as a string. Alternatively, if the stream= parameter is given a stream,
    write the file to it instead.
    """
    # If there's no stream to write to, create a StringIO to write to:
    if stream:
        xgrowf = stream
    else:
        xgrowf = StringIO()
    
    # Add some comments:
        
    for n,tile in enumerate(stxg['tiles']):
        tileedges = [str(x) for x in tile['edges']]
        for e in tileedges:
            if e != '0' and e not in [x['name'] for x in stxg['bonds']]:
                stxg['bonds'].append( { 'name': e, 'strength': 1 } )
    # Start writing the tileset
    try:
        xgrowf.write("num tile types=%d\n" % len(stxg['tiles']))
        xgrowf.write("num binding types=%d\n" % len(stxg['bonds']))
    except KeyError, e:
        raise ConvertError("Missing section {0}.".format(e.message))
  
    bondhasname = [ 'name' in x.keys() for x in stxg['bonds'] ]
    
    if all(bondhasname):
        xgrowf.write("binding type names={ %s }\n" % " ".join([x['name'] for x in stxg['bonds']]))
    elif any(bondhasname):
        raise ConvertError("Some bonds have names, but not all. All bonds must be either named or unnamed.")
    
    xgrowf.write("tile edges={\n")
    
    for n,tile in enumerate(stxg['tiles']):
        xgrowf.write("{ %s }" % " ".join([str(x) for x in tile['edges']]))
        if 'stoic' in tile:
            xgrowf.write("[%g]" % tile['stoic'])
        if 'color' in tile:
            xgrowf.write("(%s)" % tile['color'])
        xgrowf.write("   % ")
        if 'name' in tile:
            xgrowf.write(tile['name'])
        xgrowf.write(" (tile #%d)\n" % n)
    
    xgrowf.write("}\n")
    
    xgrowf.write("binding strengths={ %s }\n" % " ".join([str(x['strength']) for x in stxg['bonds']]))
    
    if 'glues' in stxg:
        for x1,x2,g in [tuple(x) for x in stxg['glues']]:
            # dereference
            if type(x1) != int:
                try:
                    x1 = [n for n,z in enumerate(stxg['bonds']) if z['name']==x1][0]+1
                except IndexError:
                    print x2
                    raise ValueError(x2)
            if type(x2) != int:
                try:
                    x2 = [n for n,z in enumerate(stxg['bonds']) if z['name']==x2][0]+1
                except IndexError:
                    print x2
                    raise ValueError(x2)
            xgrowf.write("g(%d,%d)=%g\n" % (x1,x2,float(g)))
    
    if 'xgrowargs' in stxg.keys():
        for key,val in stxg['xgrowargs'].items():
            if key in keyopts:
                xgrowf.write("%s=%s\n" % (key, str(val)))
            elif key in truefalseoptions:
                if val:
                    xgrowf.write(key+'\n')
            elif key == 'window':
                if not val:
                    xgrowf.write('-nw\n')
            elif key == 'doubletiles':
                for x1,x2 in val:
                    # dereference
                    if type(x1) != int:
                        try:
                            x1 = [n for n,z in enumerate(stxg['tiles']) if z['name']==x1][0]+1
                        except IndexError:
                            print x1
                            raise ValueError(x1)
                    if type(x2) != int:
                        try:
                            x2 = [n for n,z in enumerate(stxg['tiles']) if z['name']==x2][0]+1
                        except IndexError:
                            print x2
                            raise ValueError(x2)           
                    xgrowf.write( "doubletile=%d,%d\n" % (x1,x2) )
            elif key == 'vdoubletiles':
                for x1,x2 in val:
                    # dereference
                    if type(x1) != int:
                        try:
                            x1 = [n for n,z in enumerate(stxg['tiles']) if z['name']==x1][0]+1
                        except IndexError:
                            print x1
                            raise ValueError(x1)
                    if type(x2) != int:
                        try:
                            x2 = [n for n,z in enumerate(stxg['tiles']) if z['name']==x2][0]+1
                        except IndexError:
                            print x2
                            raise ValueError(x2)           
                    xgrowf.write( "vdoubletile=%d,%d\n" % (x1,x2) )
            elif key == 'fission':
                if val == 'off':
                    xgrowf.write("no_fission\n")
                elif val == 'on':
                    xgrowf.write("fission\n")
                elif val == 'chunk':
                    xgrowf.write("chunk_fission\n")
    xgrowf.write("% Tileset created by stxg.py version {0} on {1}\n".format(version, str(datetime.date.today())) )
    
    if stream:
        return None
    else:
        return xgrowf.getvalue()
        
        
def from_xgrow( xgst ):
    """
    Given a string of an xgrow tile file, convert it to stxg.
    
    Note that the conversion may require some work afterward.
    """
    
    ntre = re.compile(r"^\s*num tile types=(\d+)", re.MULTILINE)
    nbre = re.compile(r"^\s*num binding types=(\d+)", re.MULTILINE)
    bsre = re.compile(r"^\s*binding strengths\s*=\s*{([^}]+)}", re.MULTILINE)
    tere = re.compile(r"^\s*tile edges\s*=[\s\n]*{((?:[^{}]*{[^{}]+}[^{}]*)+)}", re.MULTILINE)
    bnre = re.compile(r"^\s*binding type names\s*=\s*{([^}]+)}", re.MULTILINE)
    glre = re.compile(r"^\s*g\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*=\s*(\S+)",re.MULTILINE)
    
    xgs = re.sub( r"%[^\n]\n", "\n", xgst )
    
    m = ntre.search(xgs)
    if not m:
        #warning.warn("Expected number of tile types.")
        pass
    ntiles = int( m.group(1) )
        
    m = nbre.search(xgs)
    if not m:
        #warning.warn("Expected number of binding types.")
        pass
    nbind = int( m.group(1) )
    
    bsstring = bsre.search(xgs).group(1)
    tstring = tere.search(xgs).group(1)
    bnmatch = bnre.search(xgs)
    if bnmatch:
        bnstring = bnmatch.group(1)
        bnames = bnstring.split()
    
    bstrengths = bsstring.split()
    tilevals = re.findall(r"{\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*}(?:\[(\S+)\])?(?:\((\S+)\))?",tstring)
    
    tiles = []
    for tv in tilevals:
        tile = {}
        tile['edges'] = []
        for e in tv[0:4]:
            try:
                e = int(e)
            except:
                pass
            tile['edges'].append(e)
        if tv[4]!='':
            tile['stoic'] = float(tv[4])
        if tv[5]!='':
            tile['color'] = tv[5]
        tiles.append(tile)
        
    if bnmatch:
        assert len(bnames) == len(bstrengths)
        bonds = [ {'name': name, 'strength': float(strength)} for name, strength in zip(bnames,bstrengths) ]
    else:
        bonds = [ {'strength': float(strength)} for strength in bstrengths ]
    
    assert len(bonds) == nbind
    assert len(tiles) == ntiles
    
    ts = {'tiles': tiles, 'bonds': bonds}
    
    glues = {}
    m = glre.findall(xgs)
    if m:
        for x1,x2,s in m:
            if bnmatch:
                glues[(bnames[int(x1)-1], bnames[int(x2)-1])] = float(s)
            else:
                glues[int(x1),int(x2)] = float(s)
        ts['glues'] = glues
    
    xgrowargs = {}
    for opt,opttype in keyvaloptions:
        m = re.findall( r"^\s*"+opt+r"\s*=\s*(\S+)", xgs, re.MULTILINE )
        if len(m)>0:
            xgrowargs[opt] = opttype(m[-1])
            
    # fission
    m = re.findall( r"^\s*(chunk_fission|fission|no_fission)", xgs, re.MULTILINE )
    if len(m)>0:
        if m[-1] == 'chunk_fission':
            xgrowargs['fission'] = 'chunk'
            
        elif m[-1] == 'fission':
            xgrowargs['fission'] = 'on'
        else:
            xgrowargs['fission'] = 'off'
    
    # doubletile
    m = re.findall( r"^\s*doubletile\s*=\s*(\d+)\s*,\s*(\d+)", xgs, re.MULTILINE )
    if m:
        xgrowargs['doubletiles'] = [ [int(x1),int(x2)] for x1,x2 in m ]
    
    if re.search(r"^\s*-nw", xgs, re.MULTILINE ):
        xgrowargs['window'] = False
        
    ts['xgrowargs'] = xgrowargs
    return ts  
    
#

# Stuff originally from yamltostxg.py. This is even less polished than everything above.

def from_yaml_endadj( ts, perfect=False, rotate=False ):
    import stickydesign as sd
    import stickydesign.energetics as en
    import numpy as np
    
    # Combine ends and tile-specified adjacents
    newtiles = []
    newends = []
    endslist = set()
    doubleends = []
    doubles = []
    vdoubleends = []
    vdoubles = []
    
    newtiles.append( { 'name': 'origami', 'edges': ['origami','origami','origami','origami'], 'stoic': 0, 'color': 'white'} )
    
    for tile in ts['seed']['adapters']:
        newtile = {}
        newtile['edges'] = [ 'origami' ] +  [ re.sub('/','_c',x) for x in tile['ends'] ] + [ 'origami' ]
        if 'name' in tile: newtile['name'] = tile['name']
        newtile['stoic'] = 0
        newtile['color'] = 'white'
        newtiles.append(newtile)    

    if rotate:
        rotatedtiles = []
        for tile in ts['tiles']:
            if tile['type'] == 'tile_daoe_3up' or tile['type'] == 'tile_daoe_5up':
                newtile = copy.copy(tile)
                newtile['name']+='_lrf'
                newtile['ends']=[tile['ends'][x] for x in (1,0,3,2)]
                rotatedtiles.append(newtile)
                newtile = copy.copy(tile)
                newtile['name']+='_udf'
                newtile['type']='tile_daoe_'+{'5up':'3up','3up':'5up'}[tile['type'][-3:]]
                newtile['ends']=[tile['ends'][x] for x in (3,2,1,0)]
                rotatedtiles.append(newtile)
                newtile = copy.copy(tile)
                newtile['name']+='_bf'
                newtile['type']='tile_daoe_'+{'5up':'3up','3up':'5up'}[tile['type'][-3:]]
                newtile['ends']=[tile['ends'][x] for x in (2,3,0,1)]
                rotatedtiles.append(newtile)
            elif tile['type'] == 'tile_daoe_doublehoriz_35up':
                newtile = copy.copy(tile)
                newtile['name']+='_lrf'
                newtile['type']='tile_daoe_doublevert_53up'
                newtile['ends']=[tile['ends'][x] for x in (2,1,0,5,4,3)]
                rotatedtiles.append(newtile)
                newtile = copy.copy(tile)
                newtile['name']+='_udf'
                newtile['type']='tile_daoe_doublevert_53up'
                newtile['ends']=[tile['ends'][x] for x in (5,4,3,2,1,0)]
                rotatedtiles.append(newtile)
                newtile = copy.copy(tile)
                newtile['name']+='_bf'
                newtile['ends']=[tile['ends'][x] for x in (3,4,5,0,1,2)]
                rotatedtiles.append(newtile)
            elif tile['type'] == 'tile_daoe_doublevert_35up':
                newtile = copy.copy(tile)
                newtile['name']+='_lrf'
                newtile['type']='tile_daoe_doublehoriz_53up'
                newtile['ends']=[tile['ends'][x] for x in (2,1,0,5,4,3)]
                rotatedtiles.append(newtile)
                newtile = copy.copy(tile)
                newtile['name']+='_udf'
                newtile['type']='tile_daoe_doublehoriz_53up'
                newtile['ends']=[tile['ends'][x] for x in (5,4,3,2,1,0)]
                rotatedtiles.append(newtile)
                newtile = copy.copy(tile)
                newtile['name']+='_bf'
                newtile['ends']=[tile['ends'][x] for x in (3,4,5,0,1,2)]
                rotatedtiles.append(newtile)

        ts['tiles'] += rotatedtiles
    
    for tile in ts['tiles']:
        if tile['type'] == 'tile_daoe_3up' or tile['type'] == 'tile_daoe_5up':
            newtile = {}
            newtile['edges'] = [ re.sub('/','_c',x) for x in tile['ends'] ]
            if 'name' in tile: newtile['name'] = tile['name']
            if 'conc' in tile: newtile['stoic'] = tile['conc']
            if 'color' in tile: newtile['color'] = tile['color']
            newtiles.append(newtile)

        if tile['type'] == 'tile_daoe_doublehoriz_35up' or tile['type'] == 'tile_daoe_doublehoriz_53up':
            newtile1 = {}
            newtile2 = {}
            newtile1['edges'] = [ re.sub('/','_c',x) for x in tile['ends'][0:1] ] \
                + [ tile['name']+'_db' ] \
                + [ re.sub('/','_c',x) for x in tile['ends'][4:] ]
            newtile2['edges'] = [ re.sub('/','_c',x) for x in tile['ends'][1:4] ] \
                + [ tile['name']+'_db' ]            
            newtile1['name'] = tile['name']+'_left'
            newtile2['name'] = tile['name']+'_right'
                        
            doubleends.append( tile['name']+'_db' )
            doubles.append( (newtile1['name'], newtile2['name']) )
            
            if 'conc' in tile: 
                newtile1['stoic'] = tile['conc']
                newtile2['stoic'] = tile['conc']
                
            if 'color' in tile: 
                newtile1['color'] = tile['color']
                newtile2['color'] = tile['color']
                
            newtiles.append(newtile1)
            newtiles.append(newtile2)
        if tile['type'] == 'tile_daoe_doublevert_35up' or tile['type'] == 'tile_daoe_doublevert_53up':
            newtile1 = {}
            newtile2 = {}
            newtile1['edges'] = [ re.sub('/','_c',x) for x in tile['ends'][0:2] ] \
                + [ tile['name']+'_db' ] \
                + [ re.sub('/','_c',x) for x in tile['ends'][5:] ]
            newtile2['edges'] = [ tile['name']+'_db' ] + [ re.sub('/','_c',x) for x in tile['ends'][2:5] ] 
            newtile1['name'] = tile['name']+'_top'
            newtile2['name'] = tile['name']+'_bottom'
                        
            vdoubleends.append( tile['name']+'_db' )
            vdoubles.append( (newtile1['name'], newtile2['name']) )
            
            if 'conc' in tile: 
                newtile1['stoic'] = tile['conc']
                newtile2['stoic'] = tile['conc']
                
            if 'color' in tile: 
                newtile1['color'] = tile['color']
                newtile2['color'] = tile['color']
                
            newtiles.append(newtile1)
            newtiles.append(newtile2)
    
    newends.append( { 'name': 'origami', 'strength': 100 } )
 
    for end in doubleends:
        newends.append( { 'name': end, 'strength': 10 } )
    for end in vdoubleends:
        newends.append( { 'name': end, 'strength': 10 } )

    gluelist = []
    if not perfect: 
        glueends = {'DT': [], 'TD': []}
        for end in ts['ends']:
            newends.append( { 'name': end['name'], 'strength': 0 } )
            newends.append( { 'name': end['name']+'_c', 'strength': 0 } )
            if (end['type'] == 'TD') or (end['type'] == 'DT'):
                glueends[end['type']].append((end['name'],end['fseq']))
                
        ef = en.energetics_santalucia(mismatchtype='max')
        
        for t in ['DT','TD']:
            names, fseqs = zip(*glueends[t])
            allnames = names + tuple( x+'_c' for x in names )
            ea = sd.endarray(fseqs, t)
            ar = sd.energy_array_uniform(ea,ef)
            for i1,n1 in enumerate(names):
                for i2,n2 in enumerate(allnames):
                    gluelist.append([n1,n2,float(ar[i1,i2])])
    else:
        if 'ends' not in ts.keys():
            ts['ends']=[]
        endsinlist = set( e['name'] for e in ts['ends'] )
        endsintiles = set()
        for tile in ts['tiles']:
            endsintiles.update( re.sub('/','',e) for e in tile['ends'] if e != 'hp')
        for end in ts['ends'] + list({'name': e} for e in endsintiles):
            newends.append( { 'name': end['name'], 'strength': 0 } )
            newends.append( { 'name': end['name']+'_c', 'strength': 0 } )
            gluelist.append([end['name'],end['name']+'_c',1.0]) 
            

    
    newends.append( {'name': 'hp', 'strength': 0} )

    xga = {}
    xga['doubletiles'] = [ list(x) for x in doubles ]
    xga['vdoubletiles'] = [ list(x) for x in vdoubles ]
    xga.update( ts['xgrow_options'] )
    xga.update( ts['xgrow_options'] )
     
        
    sts = { 'tiles': newtiles, 'bonds': newends, 'xgrowargs': xga, 'glues': gluelist }
    
    return sts

def from_yaml_tileadj( ts ):
    import stickydesign as sd
    import numpy as np

    # Combine ends and tile-specified adjacents
    newtiles = []
    newends = []
    endslist = set()
    doubleends = []
    doubles = []
    
    newtiles.append( { 'name': 'origami', 'edges': ['origami','origami','origami','origami'], 'stoic': 1e-9, 'color': 'white'} )
    
    for tile in ts['seed']['adapters']:
        newtile = {}
        newtile['edges'] = [ 'origami' ] +  [ re.sub('/','_c',x)+'_'+y for x,y in zip(tile['ends'],tile['adjs']) ] + [ 'origami' ]
        endslist.update( set(zip( newtile['edges'][1:3], tile['ends'], tile['adjs'] )) )
        if 'name' in tile: newtile['name'] = tile['name']
        newtile['stoic'] = 1e-9
        newtile['color'] = 'white'
        newtiles.append(newtile)    

    for tile in ts['tiles']:
        if tile['type'] == '3up' or tile['type'] == '5up':
            newtile = {}
            newtile['edges'] = [ re.sub('/','_c',x)+'_'+y for x,y in zip(tile['ends'],tile['adjs']) ]
            endslist.update( set(zip( newtile['edges'], tile['ends'], tile['adjs'] )) )
            if 'name' in tile: newtile['name'] = tile['name']
            if 'conc' in tile: newtile['stoic'] = tile['conc']
            if 'color' in tile: newtile['color'] = tile['color']
            newtiles.append(newtile)
        if tile['type'] == '3up5up' or tile['type'] == '5up3up':
            newtile1 = {}
            newtile2 = {}
            newtile1['edges'] = [ re.sub('/','_c',x)+'_'+y for x,y in zip(tile['ends'][0:1],tile['adjs'][0:1]) ] \
                + [ tile['name']+'_db' ] \
                + [ re.sub('/','_c',x)+'_'+y for x,y in zip(tile['ends'][4:],tile['adjs'][4:]) ]
            newtile2['edges'] = [ re.sub('/','_c',x)+'_'+y for x,y in zip(tile['ends'][1:4],tile['adjs'][1:4]) ] \
                + [ tile['name']+'_db' ]
            endslist.update( set(zip( [ re.sub('/','_c',x)+'_'+y for x,y in zip(tile['ends'],tile['adjs']) ], tile['ends'], tile['adjs'] )) )
            
            newtile1['name'] = tile['name']+'_left'
            newtile2['name'] = tile['name']+'_right'
                        
            doubleends.append( tile['name']+'_db' )
            doubles.append( (newtile1['name'], newtile2['name']) )
            
            if 'conc' in tile: 
                newtile1['stoic'] = tile['conc']
                newtile2['stoic'] = tile['conc']
                
            if 'color' in tile: 
                newtile1['color'] = tile['color']
                newtile2['color'] = tile['color']
                
            newtiles.append(newtile1)
            newtiles.append(newtile2)
    
    newends.append( { 'name': 'origami', 'strength': 100 } )
    
    for end in doubleends:
        newends.append( { 'name': end, 'strength': 100 } )
        
    for end in endslist:
        newends.append( { 'name': end[0], 'strength': 0 } )
    
    xga = {}
    xga['doubletiles'] = [ list(x) for x in doubles ]
    xga.update( ts['xgrow_options'] )
    
    
    # Now finally, deal with glues... in a perfect world first?
    # gluelist = []
    # for end1 in endslist:
    #     if end1[1][-1] == '/': continue
    #     for end2 in ( x for x in endslist if x[1]==end1[1]+'/' ):
    #         gluelist.append( [end1[0],end2[0],1] )
    #         
    ef = sd.energyfuncs_santalucia(mismatchtype='max')
    # 
    # gluelist = []
    # for end1 in endslist:
    #     if end1[1][-1] != '/':
    #         ed1 = [ x for x in ts['ends'] if x['name'] == end1[1] ][0]
    #     else:
    #         edd1 = [ x for x in ts['ends'] if x['name'] == end1[1][:-1] ][0]
    #         ed1 = { 'type': edd1['type'], 'seq': revcomp(edd1['seq'])}
    #     if (ed1['type'] == 'fake') or (ed1['type'] == 'fakeDT') or (ed1['type'] == 'fakeTD'):
    #         continue
    #     elif (ed1['type'] == 'DT'):
    #         e1 = sd.endarray([(end1[2]+ed1['seq']).lower()],'DT')
    #     elif (ed1['type'] == 'TD'):
    #         e1 = sd.endarray([(ed1['seq']+end1[2]).lower()],'TD')
    #     for end2 in endslist:
    #         if end2[1][-1] != '/':
    #             ed2 = [ x for x in ts['ends'] if x['name'] == end2[1] ][0]
    #         else:
    #             edd2 = [ x for x in ts['ends'] if x['name'] == end2[1][:-1] ][0]
    #             ed2 = { 'type': edd2['type'], 'seq': revcomp(edd2['seq'])}
    #         if ed1['type'] != ed2['type']:
    #             continue
    #         elif (ed1['type'] == 'DT'):
    #             e2 = sd.endarray([(end2[2]+ed2['seq']).lower()],'DT')
    #         elif (ed1['type'] == 'TD'):
    #             e2 = sd.endarray([(ed2['seq']+end2[2]).lower()],'TD')
    #         gluelist.append( [end1[0], end2[0], float(ef.uniform(e1,e2)[0])])
    dtl1 = []
    dtl2 = []
    tdl1 = []
    tdl2 = []
    dtn1 = []
    dtn2 = []
    tdn1 = []
    tdn2 = []
    for end1 in endslist:
        if end1[1][-1] != '/':
            ed1 = [ x for x in ts['ends'] if x['name'] == end1[1] ][0]
        else:
            edd1 = [ x for x in ts['ends'] if x['name'] == end1[1][:-1] ][0]
            ed1 = { 'type': edd1['type'], 'seq': revcomp(edd1['seq'])}
        if (ed1['type'] == 'fake') or (ed1['type'] == 'fakeDT') or (ed1['type'] == 'fakeTD'):
            continue
        for end2 in endslist:
            if end2[1][-1] != '/':
                ed2 = [ x for x in ts['ends'] if x['name'] == end2[1] ][0]
            else:
                edd2 = [ x for x in ts['ends'] if x['name'] == end2[1][:-1] ][0]
                ed2 = { 'type': edd2['type'], 'seq': revcomp(edd2['seq'])}
            if ed1['type'] != ed2['type']:
                continue
            elif (ed1['type'] == 'DT'):
                dtl1.append( (end1[2]+ed1['seq']).lower() )
                dtl2.append( (end2[2]+ed2['seq']).lower() )
                dtn1.append( end1[0] )
                dtn2.append( end2[0] )
            elif (ed1['type'] == 'TD'):
                tdl1.append( (ed1['seq']+end1[2]).lower() )
                tdl2.append( (ed2['seq']+end2[2]).lower() )
                tdn1.append( end1[0] )
                tdn2.append( end2[0] )
    
    dta1 = sd.endarray( dtl1, 'DT' )
    dta2 = sd.endarray( dtl2, 'DT' )
    tda1 = sd.endarray( tdl1, 'TD' )
    tda2 = sd.endarray( tdl2, 'TD' )
    
    dtg = ef.uniform(dta1, dta2)
    tdg = ef.uniform(tda1, tda2)
    
    dtgl = set([ tuple(sorted([x,y]) + [float(z)]) for x,y,z in zip(dtn1, dtn2, dtg) ])
    tdgl = set([ tuple(sorted([x,y]) + [float(z)]) for x,y,z in zip(tdn1, tdn2, tdg) ])
    
    gluelist = [ list(x) for x in dtgl.union(tdgl) ]
        
    sts = { 'tiles': newtiles, 'bonds': newends, 'xgrowargs': xga, 'glues': gluelist }
    
    return sts

def revcomp(seq):
    import stickydesign as sd
    import numpy as np
    seq = seq.lower()
    return "".join( reversed( [sd.wc[nt] for nt in seq] ) )
    
