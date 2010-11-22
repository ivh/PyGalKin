"""
db.py

"""

from PyGalKin import *

DBNAME='/data/cigale/cigale.db'
exists=path.exists

def setupdb(dbname=DBNAME):
    if not exists(dbname): return "Database not there!"
    connection=sqlite.connect(dbname)
    cursor=connection.cursor()
    return connection,cursor

def getg(curs,cols,where='gid NOTNULL',table='galax'):
    curs.execute("SELECT %s FROM %s WHERE %s"%(cols,table,where))
    return map(list,curs.fetchall())

def geto(curs,cols,where='oid NOTNULL',table='obspar'):
    curs.execute("SELECT %s FROM %s WHERE %s"%(cols,table,where))
    return N.transpose(map(list,curs.fetchall()))

def getd(curs,cols,where='did NOTNULL',table='data'):
    curs.execute("SELECT %s FROM %s WHERE %s"%(cols,table,where))
    return N.transpose(map(list,curs.fetchall()))

def getSDSSids(curs):
    tol=0.00833
    ids=getg(curs,'gid')[0]
    ra,dec=getg(curs,'ra,dec')
    names=getg(curs,'name')[0]
    for i,name in enumerate(names):
        dec0,dec1=dec[i]-tol,dec[i]+tol
        ra0,ra1=ra[i]-tol,ra[i]+tol
        if True:
        #if "UM 501" in name:
            print i,name,ra[i],dec[i]
            print sqlcl.sqlcl("SELECT objid,specobjid,ra,dec,z from SpecPhotoAll WHERE (ra BETWEEN %s AND %s) AND (dec BETWEEN %s AND %s) AND (specobjid !=0) AND (primtarget & 0x40 > 0)"%(ra0,ra1,dec0,dec1))
