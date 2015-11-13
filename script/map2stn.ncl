; Written by Wang Wenshan 2014-07-01 Tuesday 09:55:29
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;*************************************************************
; First step to process AIRS data
;   - var (hdf --> netcdf)
;   - stn
;       set buffer: min of all map doesn't mean close to stn
;       read 100 min together: stn at edges
;   - pressure reverse
;
; To-do (July 6, 2014; done by July 7, 2014):
;   - use only land
;   - combine nomiss.ncl (not the whole script but chose non-missing pfl)
;
; Add: Nov. 6, 2014
;   - more var: totH2OStd
;   - use both QC and error
;
; Resource: map_lat-stn.res
;************************************************************
begin
;- tuning
nmiss=3

;- para
yr="2010"
;yr="2010.07.03"       ; swiss only now
fname = "map2stn"
buf = 25.*sqrt(2)       ; buffer around one gc-net stn in km 
;- the furthest one is just at one corner of the AIRS grid box

nx = 45
ny = 30
npt = 28        ; plev of tmp
npq = 15        ; plev of H2O

;- AIRS files
drc = "/data/wenshanw/airs/"
drc_in = drc+"hdf/"
drc_out = drc+"stn/"
;drc_out = "/data/wenshanw/crm/xpt/run/"
fns = systemfunc("ls "+drc_in+"AIRS."+yr+".*.hdf")
tmid = toint(systemfunc("ls "+drc_in+"AIRS."+yr+"*.hdf | cut -d '.' -f 5"))
date = systemfunc("ls "+drc_in+"AIRS."+yr+"*.hdf | cut -d '.' -f 2-5")
nfn = dimsizes(fns)
;print(""+fns)
;print(""+tmid)
;print(""+date)
;exit
ifn = 0       ; total number of file processed

;- GC-Net stations
stns = systemfunc("cat ~/cloud/doc/stn/stn14_crd_gc-net.txt | awk '{print $1}'")
latstn = tofloat(systemfunc("cat ~/cloud/doc/stn/stn14_crd_gc-net.txt | awk '{print $2}'"))
lonstn = tofloat(systemfunc("cat ~/cloud/doc/stn/stn14_crd_gc-net.txt | awk '{print $3}'"))
;stns = systemfunc("cat ~/cloud/doc/stn/stn_humb.txt | awk '{print $1}'")
;latstn = tofloat(systemfunc("cat ~/cloud/doc/stn/stn_humb.txt | awk '{print $2}'"))
;lonstn = tofloat(systemfunc("cat ~/cloud/doc/stn/stn_humb.txt | awk '{print $3}'"))
nstn = dimsizes(stns)
latstn3d = conform_dims((/nstn,nx,ny/),latstn,(/0/))
lonstn3d = conform_dims((/nstn,nx,ny/),lonstn,(/0/))
latstn3d@units = "degrees_north"
lonstn3d@units = "degrees_east"
;print(""+nstn)
;print(lonstn+" "+latstn)
;exit

;- var
;- renew every stn every time point
;v3d = (/"TAirStd","H2OMMRLevStd"/)
v3d = (/"TAirStd","H2OMMRLevStd","O3VMRLevStd"/)
v2d = (/"TSurfAir","H2OMMRSurf","PSurfStd","nSurfStd","totH2OStd","Time","solzen","solazi","landFrac","CldFrcTot"/)
;v2d = (/"PSurfStd","Time","solzen","solazi","landFrac","CldFrcTot"/)
;v2d = (/"PSurfStd","Time","solzen","solazi","CldFrcTot"/)
v3dqc = (/"TAirStd_QC","H2OMMRLevStd_QC","O3VMRLevStd_QC"/)
v2dqc = (/"TSurfAir_QC","H2OMMRSurf_QC","PSurfStd_QC","totH2OStd_QC","CldFrcTot_QC"/)
;v2dqc = (/"PSurfStd_QC"/)
v3derr = (/"TAirStdErr","H2OMMRLevStdErr","O3VMRLevStdErr"/)
v2derr = (/"TSurfAirErr","H2OMMRSurfErr","totH2OStdErr","landFrac_err"/)
v3du = (/"K","gm/kg","kg/kg"/)
v2du = (/"K","gm/kg","hPa","1","kg/m^2","seconds since 1993-1-1 00:00:00","degree","degree","1","1"/)
;v2du = (/"hPa","seconds since 1993-1-1 00:00:00","degree","degree","1","1"/)
;v2du = (/"hPa","seconds since 1993-1-1 00:00:00","degree","degree","1"/)
vqcu = "0,1,2"        ; qc units
v3derru = v3du
v2derru = (/"K","gm/kg","kg/m^2","1"/)
nv3d = dimsizes(v3d)
nv2d = dimsizes(v2d)
nv3dqc = dimsizes(v3dqc)
nv2dqc = dimsizes(v2dqc)
nv3derr = dimsizes(v3derr)
nv2derr = dimsizes(v2derr)

;- renew every time
;dst_min = new((/nstn/), double)        ; min distance for each stn
;dst_min_ind = new((/nstn,3/),integer) ; ind of min distance
dst_fn = new((/nstn/),integer)         ; which file has the min
dst_std = new((/nstn/),double)        ; min of all time
dst_std = buf 
dst_min_ind = new((/nstn,2/),integer)     ; min index
xid = new((/nstn/),integer)        ; min index of all time
yid = new((/nstn/),integer)

do while (ifn .lt. nfn)       ; loop: files
;- read 100 min file together 
;- deal with 100 min
;   update ifn
  print("Begin! ifn = "+ifn)
;  print("fn = "+fns(ifn))
;- map --> stn
  fin = addfile(fns(ifn),"r")
  lat = fin->Latitude
  dim_lat = dimsizes(lat)
;- if granule is not 45*30
  if (.not. (dim_lat(0).eq.nx.and.dim_lat(1).eq.ny)) then
    ;- clearup
    print("Out of Shp: "+fns(ifn))
    print(ifn+" "+dim_lat(0)+" "+dim_lat(1))
    delete([/lat,dim_lat/])
    continue
  end if

;- QC and land
  tqc = fin->TAirStd_QC
  qqc = fin->H2OMMRLevStd_QC
  tsqc = fin->TSurfAir_QC
  psqc = fin->PSurfStd_QC

;- land fration error: drop (Nov. 6, 2014)
  ld = fin->landFrac
;  lderr = fin->landFrac_err
;  printMinMax(ld,True)
;  print(num(.not.ismissing(ld)))
;  ld = where(lderr.lt.0.or.lderr.gt.ld*0.5,ld@_FillValue, ld)
;  printMinMax(ld,True)
;  print(num(.not.ismissing(ld)))
;  exit
;  printVarSummary(tqc)
  lat = where(dim_min(tqc).eq.2, -89., lat)        ; set the bad points to North Pole
  lat = where(dim_min(qqc).eq.2, -89., lat)        ; set the bad points to North Pole
  lat = where(tsqc.eq.2, -89., lat)        ; set the bad points to North Pole
  lat = where(psqc.eq.2, -89., lat)        ; set the bad points to North Pole
  lat = where(ld.lt.0.75, -89., lat)       ; if not land points, set to North Pole
;  print("Missing Temp Profile = "+num(dim_min(tqc).eq.2))
;  print("Missing H2O Profile = "+num(dim_min(qqc).eq.2))
;  print("Ocean Grid = "+num(ld .lt. 0.75))

;- error
  ts = fin->TSurfAir
  tserr = fin->TSurfAirErr
  lat = where(tserr.lt.0 .or. tserr.gt.ts*0.5, -89., lat)
  qs = fin->H2OMMRSurf
  qserr = fin->H2OMMRSurfErr
  lat = where(qserr.lt.0 .or. qserr.gt.qs*0.5, -89., lat)
  delete([/ts,tserr,qs,qserr/])

  t = fin->TAirStd
  terr = fin->TAirStdErr
  t = where(terr.lt.0 .or. terr.gt.t*0.5, t@_FillValue, t)
  q = fin->H2OMMRLevStd
  qerr = fin->H2OMMRLevStdErr
  q = where(qerr.lt.0 .or. qerr.gt.q*0.5, q@_FillValue, q)
  delete([/terr,qerr/])

;- actual miss
  nt = dim_num(ismissing(t(:,:,4:)))
  lat = where(nt.gt.nmiss,-89.,lat)         ; lower layer cannot miss more than nmiss 
  nq = dim_num(ismissing(q(:,:,4:)))
;  print(q&H2OPressureLev_L2_Standard_atmospheric_surface_product)
;  printMinMax(nq,True)
;  printVarSummary(q)
;  exit
  lat = where(nt.gt.nmiss,-89.,lat)         ; lower layer cannot miss more than nmiss 
  delete([/tqc,qqc,tsqc,psqc,ld,t,nt,q,nq/])

  lon = fin->Longitude
  lat3d = conform_dims((/nstn,nx,ny/),lat,(/1,2/)) 
  lon3d = conform_dims((/nstn,nx,ny/),lon,(/1,2/)) 
  lat3d@units = "degrees_north"
  lon3d@units = "degrees_east"
;- Distance Method 1: sqrt(dx^2+dy^2)
;- Distance Method 2: great circle distance
  dst = gc_latlon(lat3d,lon3d,latstn3d,lonstn3d,2,-4)
  dst_reshp = reshape(dst,(/nstn,nx*ny/))
  dst_min = dim_min(dst_reshp)        ; nstn, min distance 
  dst_min_3d = conform_dims((/nstn,nx,ny/),dst_min,(/0/))
  dst_min_ind_all = ind_resolve(ind(ndtooned(dst).eq.ndtooned(dst_min_3d)),(/nstn,nx,ny/))
  dst_min_ind(dst_min_ind_all(:,0),0) = dst_min_ind_all(:,1)
  dst_min_ind(dst_min_ind_all(:,0),1) = dst_min_ind_all(:,2)

;  print(dst_min_ind(:,0)+" "+dst_min_ind(:,1))
;  do j = 0, 12
;    print(dst(j,dst_min_ind(j,0),dst_min_ind(j,1))+" "+dst_min(j))
;    print(dst(j,xid(j),yid(j))+" "+dst_min(j))
;  end do

;- which one to keep
  xid = where(dst_min.lt.dst_std,dst_min_ind(:,0),xid)      ; min dst ind
  yid = where(dst_min.lt.dst_std,dst_min_ind(:,1),yid)
  dst_fn = where(dst_min.lt.dst_std,ifn,dst_fn)
  dst_std = where(dst_min.lt.dst_std,dst_min,dst_std)       ; min dst
;- check crd
;  print(xid+" "+yid+" "+dst_fn+" "+dst_std)
;  do jstn = 0, nstn - 1
;    if (.not.ismissing(xid(jstn)).and.dst_fn(jstn).eq.ifn) then
;      print(stns(jstn)+"("+latstn(jstn)+","+lonstn(jstn)+") --> ("+lat(xid(jstn),yid(jstn))+","+lon(xid(jstn),yid(jstn))+"): "+dst_std(jstn)+" km; ifn = "+dst_fn(jstn))
;    end if
;  end do

;- clearup every file
  delete([/lat,lon,lat3d,lon3d/])
  delete([/dst,dst_reshp,dst_min,dst_min_3d,dst_min_ind_all/])

;- whether within 100 min
  if (ifn .lt. nfn-1) then        ; if the last file; go to output directly
    if (tmid(ifn+1).eq.tmid(ifn)+1) then        ; next time id is consecutive
      ifn = ifn + 1
  ;    print("Read next file in 100 min; ifn = "+ifn)
      continue
    end if
  end if

;  print("100 min done; ifn = "+ifn)
;- assign data
  do jstn = 0, nstn - 1        ; loop: stn
    if (dst_std(jstn).ge.buf) then        ; if no airs data found for this station
;      print("No data for "+stns(jstn))
      continue
    end if

;    if (stns(jstn).ne."swiss") then
;      continue
;    end if

;- v3d = (/"TAirStd","H2OMMRLevStd"/)
;- v2d = (/"PSurfStd","Time","solzen","solazi","CldFrcTot"/)

;    print(stns(jstn)+": "+dst_std(jstn)+" km; ifn = "+dst_fn(jstn))
    fin = addfile(fns(dst_fn(jstn)),"r")
    fo = drc_out+stns(jstn)+"."+date(dst_fn(jstn))+".nc"
    system("/bin/rm -f "+fo)
    fout = addfile(fo,"c")
    do kvar = 0, nv3d - 1        ; loop: 3D var 
      dt = fin->$v3d(kvar)$(xid(jstn),yid(jstn),::-1)       ; reverse plev !!!
      dt@units = v3du(kvar)
;      dt@_FillValue = -9999
      fout->$v3d(kvar)$ = dt
;      print(v3d(kvar)+" = "+dt+" "+dt@units)
      delete(dt)
    end do        ; loop done: 3D var
    do kvar = 0, nv2d - 1        ; loop: 2D var 
      dt = fin->$v2d(kvar)$(xid(jstn),yid(jstn))
      dt@units = v2du(kvar)
      fout->$v2d(kvar)$ = dt
;      print(v2d(kvar)+" = "+dt+" "+dt@units)
      delete(dt)
    end do        ; loop done: 2D var
    do kvar = 0, nv3dqc - 1        ; loop: 3D var QC 
      tt = fin->$v3dqc(kvar)$(xid(jstn),yid(jstn),::-1)       ; reverse plev !!!
      dt = toint(tt)
      copy_VarMeta(tt,dt)
      dt@units = vqcu
      fout->$v3dqc(kvar)$ = dt
;      print(v3dqc(kvar)+" = "+dt)
      delete([/tt,dt/])
    end do        ; loop done: 3D var QC
    do kvar = 0, nv2dqc - 1        ; loop: 2D var QC
      tt = fin->$v2dqc(kvar)$(xid(jstn),yid(jstn))
      dt = toint(tt)
      copy_VarMeta(tt,dt)
      dt@units = vqcu
      fout->$v2dqc(kvar)$ = dt
;      print(v2dqc(kvar)+" = "+dt)
      delete([/tt,dt/])
    end do        ; loop done: 2D var QC
;- error
    do kvar = 0, nv3derr - 1        ; loop: 3D var QC 
      tt = fin->$v3derr(kvar)$(xid(jstn),yid(jstn),::-1)        ; reverse pressure !!!
      dt = toint(tt)
      copy_VarMeta(tt,dt)
      dt@units = v3derru(kvar)
      fout->$v3derr(kvar)$ = dt
;      print(v3dqc(kvar)+" = "+dt)
      delete([/tt,dt/])
    end do        ; loop done: 3D var QC
    do kvar = 0, nv2derr - 1        ; loop: 2D var QC
      tt = fin->$v2derr(kvar)$(xid(jstn),yid(jstn))
      dt = toint(tt)
      copy_VarMeta(tt,dt)
      dt@units = v2derru(kvar)
      fout->$v2derr(kvar)$ = dt
;      print(v2dqc(kvar)+" = "+dt)
      delete([/tt,dt/])
    end do        ; loop done: 2D var QC

  end do        ; loop done: stn, output data

;- clearup every 100 min
  xid= xid@_FillValue
  yid= yid@_FillValue
  dst_std= buf 
  dst_fn = dst_fn@_FillValue

  ifn = ifn + 1

;  if (ifn .gt. 7) then      ; TEMP: try file 0-7
;    break
;  end if
end do        ; loop done: files

print("done!")
end
