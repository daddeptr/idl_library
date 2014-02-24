pro merge_xfaster_newdat_files, file1, file2, l=l, bin=bin, xfdir=xfdir, old=old, outfile=outfile

    if (not keyword_set(outfile)) then outfile = xfdir+'merge.dat'
    if (not keyword_set(xfdir)) then xfdir = '/global/scratch2/sd/dpietrob/Software/XFaster/'
    if (not keyword_set(bin)) and (not keyword_set(l)) then stop, 'ERROR: Neither l nor bin specified.'

;## NB consistent with my definition of .newdat
    if keyword_set(old) then readcol,file1,nbtt1,nbee1,nbbb1,nbtb1,nbte1,nbeb1, format='f,f,f,f,f,f',skipline=2,numline=1 else $
      readcol,file1,nbtt1,nbee1,nbbb1,nbte1,nbtb1,nbeb1, format='f,f,f,f,f,f',skipline=2,numline=1
        
    if keyword_set(old) then readcol,file2,nbtt2,nbee2,nbbb2,nbtb2,nbte2,nbeb2, format='f,f,f,f,f,f',skipline=2,numline=1 else $
      readcol,file2,nbtt2,nbee2,nbbb2,nbte2,nbtb2,nbeb2, format='f,f,f,f,f,f',skipline=2,numline=1
        
    print, ' TT # bin = ', nbtt1, nbtt2
    print, ' EE # bin = ', nbee1, nbee2
    print, ' BB # bin = ', nbbb1, nbbb2
    print, ' TE # bin = ', nbte1, nbte2
    print, ' TB # bin = ', nbtb1, nbtb2
    print, ' EB # bin = ', nbeb1, nbeb2
    
    openw,1,outfile+'.tmp'
; ------- TT
    printf,1,'TT'
    readcol,file1,cl1,cler1,lmn1,lmx1,format='x,f,f,x,x,f,f', skipline=7, numline=nbtt1[0]
    readcol,file2,cl2,cler2,lmn2,lmx2,format='x,f,f,x,x,f,f', skipline=7, numline=nbtt2[0]
    for i=0l,nbtt1[0]-1 do begin
        if lmn1[i] lt l then printf,1,cl1[i],cler1[i],lmn1[i],lmx1[i]
    endfor
    for i=0l,nbtt2[0]-1 do begin
        if lmn2[i] ge l then printf,1,cl2[i],cler2[i],lmn2[i],lmx2[i]
    endfor

; ------- EE
    printf,1,'EE'
    readcol,file1,cl1,cler1,lmn1,lmx1,format='x,f,f,x,x,f,f', skipline=7+1+2*nbtt1[0], numline=nbee1[0]
    readcol,file2,cl2,cler2,lmn2,lmx2,format='x,f,f,x,x,f,f', skipline=7+1+2*nbtt2[0], numline=nbee2[0]
    for i=0l,nbee1[0]-1 do begin
        if lmn1[i] lt l then printf,1,cl1[i],cler1[i],lmn1[i],lmx1[i]
    endfor
    for i=0l,nbee2[0]-1 do begin
        if lmn2[i] ge l then printf,1,cl2[i],cler2[i],lmn2[i],lmx2[i]
    endfor

; ------- BB
    printf,1,'BB'
    readcol,file1,cl1,cler1,lmn1,lmx1,format='x,f,f,x,x,f,f', skipline=7+1+2*nbtt1[0]+2*nbee1[0]+1, numline=nbbb1[0]
    readcol,file2,cl2,cler2,lmn2,lmx2,format='x,f,f,x,x,f,f', skipline=7+1+2*nbtt2[0]+2*nbee2[0]+1, numline=nbbb2[0]
    for i=0l,nbbb1[0]-1 do begin
        if lmn1[i] lt l then printf,1,cl1[i],cler1[i],lmn1[i],lmx1[i]
    endfor
    for i=0l,nbbb2[0]-1 do begin
        if lmn2[i] ge l then printf,1,cl2[i],cler2[i],lmn2[i],lmx2[i]
    endfor

; ------- TE
    printf,1,'TE'
    readcol,file1,cl1,cler1,lmn1,lmx1,format='x,f,f,x,x,f,f', skipline=7+1+2*nbtt1[0]+2*nbee1[0]+1+2*nbbb1[0]+1, numline=nbte1[0]
    readcol,file2,cl2,cler2,lmn2,lmx2,format='x,f,f,x,x,f,f', skipline=7+1+2*nbtt2[0]+2*nbee2[0]+1+2*nbbb2[0]+1, numline=nbte2[0]
    for i=0l,nbte1[0]-1 do begin
        if lmn1[i] lt l then printf,1,cl1[i],cler1[i],lmn1[i],lmx1[i]
    endfor
    for i=0l,nbte2[0]-1 do begin
        if lmn2[i] ge l then printf,1,cl2[i],cler2[i],lmn2[i],lmx2[i]
    endfor

; ------- TB
    printf,1,'TB'
    readcol,file1,cl1,cler1,lmn1,lmx1,format='x,f,f,x,x,f,f', skipline=7+1+2*nbtt1[0]+2*nbee1[0]+1+2*nbbb1[0]+1+2*nbte1[0]+1, numline=nbtb1[0]
    readcol,file2,cl2,cler2,lmn2,lmx2,format='x,f,f,x,x,f,f', skipline=7+1+2*nbtt2[0]+2*nbee2[0]+1+2*nbbb2[0]+1+2*nbte2[0]+1, numline=nbtb2[0]
    for i=0l,nbtb1[0]-1 do begin
        if lmn1[i] lt l then printf,1,cl1[i],cler1[i],lmn1[i],lmx1[i]
    endfor
    for i=0l,nbtb2[0]-1 do begin
        if lmn2[i] ge l then printf,1,cl2[i],cler2[i],lmn2[i],lmx2[i]
    endfor

; ------- EB
    printf,1,'EB'
    readcol,file1,cl1,cler1,lmn1,lmx1,format='x,f,f,x,x,f,f', skipline=7+1+2*nbtt1[0]+2*nbee1[0]+1+2*nbbb1[0]+1+2*nbte1[0]+1+2*nbtb1[0]+1, numline=nbeb1[0]
    readcol,file2,cl2,cler2,lmn2,lmx2,format='x,f,f,x,x,f,f', skipline=7+1+2*nbtt2[0]+2*nbee2[0]+1+2*nbbb2[0]+1+2*nbte2[0]+1+2*nbtb2[0]+1, numline=nbeb2[0]
    for i=0l,nbeb1[0]-1 do begin
        if lmn1[i] lt l then printf,1,cl1[i],cler1[i],lmn1[i],lmx1[i]
    endfor
    for i=0l,nbeb2[0]-1 do begin
        if lmn2[i] ge l then printf,1,cl2[i],cler2[i],lmn2[i],lmx2[i]
    endfor
    close,1

spawn, 'head '+outfile+'.tmp'
spawn, 'tail '+outfile+'.tmp'

readcol, outfile+'.tmp', t
openw,1,'head.tmp'
printf,1,outfile
printf,1,'TT EE BB TE TB EB'
cnt = 0

end
