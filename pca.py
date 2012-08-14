#Python script to take data from raw format or phased output and run smartpca
#Set up test case from raw data to phased to plot, just for Angus and Brahman
#BE FUCKING CAREFUL ABOUT SNP ID  NUMBERS


import os
import glob
import copy
import operator
import pickle
import random
from time import sleep


def make_dicts():
##Locdict will take a SNP Name and output a umd3 SNP number
    locdict={} 
    fil=open('''SNP_Map.txt''', 'r')
    for lin in fil:
        row=lin.split()
        locdict[row[1][1:-1]]=row[0]
    fil.close()
    ibmclocdict={} ##takes an ibmc number (AS IN RAW DATA FILES)  and outputs SNP name
    fil=open("ibmc_markers_041709.csv", 'r')
    for lin in fil:
        row=lin.split(',')
        num=row[1]
        nam=row[9]
        ibmclocdict[num]=nam
    fil.close()
    forw_dict={}
    for item in ibmclocdict.keys():
        try:
            forw_dict[item]=locdict[ibmclocdict[item]]
        except: pass
    rev_dict={} #takes SNP_Map #, and turns it back into icbm_locdict #
    for item in forw_dict.keys():
        rev_dict[forw_dict[item]]=item
    #this section makes a dictionary of all SNPs on each chromosome
    chrm_dict={}     #it is isn SNP map numbers.
    chrm_n=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29","30"]
    for item in chrm_n:
        chrm_dict[item]=[]
    nam_dict={}#WTF does namdict do?
    locations=[]
    markers=[]
    chrms=[]
    fil=open("SNP_Map.txt", 'r')
    fi=fil.readlines()
    fi=fi[1:]
    for lin in fi:
          row=lin.split('\t')
          num=int(row[0])
          nam=row[1].strip('"')
          nam_dict[nam]=num
          chrm=row[2]
          if chrm =='"X"':
                chrm='30'
          pos=int(row[3])
          if chrm in chrm_n and pos != 0:
                locs=[num, nam, chrm, pos]
                locations.append(locs)
                chrms.append(chrm)
                chrm_dict[chrm].append(locs)
          else:
                pass
    for item in chrm_dict.keys():
      chrm_dict[item].sort(key=operator.itemgetter(3))
        ##'''if num in colsort: ##this is key- excludes SNPS that are not in the struct data file'''
    trans={}#translates from Raw to number of alleles of 2
    trans['1']='0'
    trans['2']='2'
    trans['3']='1'
    trans['10']='10'
    ind_info=open('ids.csv','r').readlines()
    inddict={}
    for lin in ind_info:
        lii=lin.split(',')
        inddict[lii[0]]=lii[1:]
    ABdict={}
    ABdict['1','1']='1'
    ABdict['2','2']='2'
    ABdict['1','2']='3'
    ABdict['2','1']='3'
    return (ibmclocdict, locdict, nam_dict,trans,chrm_dict,chrm_n,forw_dict,rev_dict,inddict,ABdict)


def make3k_revdict():
    (ibmclocdict, locdict, nam_dict,trans,chrm_dict,chrm_n,forw_dict,rev_dict,inddict,ABdict)=make_dicts()
##Locdict will take a SNP_Map3k Name and output a umd3 SNP number
    locdict3k={}
    fil=open('''SNP_Map3k.txt''', 'r')
    for lin in fil:
        row=lin.split()
        locdict3k[row[1]]=row[0]
    fil.close()
    forw_dict3k={}#takes ibmc # and gives SNP_Map3k number
    for item in ibmclocdict.keys():
        try:
            forw_dict3k[item]=locdict3k[ibmclocdict[item]]
        except: pass
    rev_dict3k={} #takes SNP_Map #, and turns it back into icbm_locdict #
    for item in forw_dict3k.keys():
        rev_dict3k[forw_dict3k[item]]=item
    return rev_dict3k


def raw_to_struct(finame):
    (ibmclocdict, locdict, nam_dict,trans,chrm_dict,chrm_n,forw_dict,rev_dict,inddict,ABdict)=make_dicts()
    print finame
    fi=open(finame)
    rawdat=[]
    i=1
    for lin in fi:
            stri=str(str(lin)[1:-2])
            stri2=(str(stri)).split(',')
            try:
                #renames snpnum from ibmc numbers to SNP_Map numbers in rawdat
                snpnum=nam_dict[ibmclocdict[stri2[0].strip()]] 
                lst=[snpnum,stri2[1].strip(),stri2[2].strip()]
                rawdat.append(lst)                                                    
            except:
                i=i+1
    print "Skipped lines = %i" %i
    fi.close()
    print 'data read'
    nudat=rawdat
    #now nudat holds all the raw data recoded as [SNP_Map snpnumber,individual id number, genotype]
    snp_dict={}
    indiv_n=set()
    for item in nudat:#dictionarizes nudat
        indiv_n.add(item[1])#makes a lis of inidvdual ids
        snp_dict[str(item[0]),item[1]]=item[2]#tuple key of (SNP_Map #, ID #) calls genotype
    
    indiv_n=list(indiv_n)
    indiv_n.sort()
    
    
    for item in chrm_n: #goes through each chromosome SNP dictionary, sorts it, and then pulls the genotypes for that SNP.
        print nm
        #here globals()[nm] calls the list called chrm1 then sorts it
        #return locals()[nm]
        chrm_dict[chrm_n]=sorted(chrm_dict[chrm_n], key=operator.itemgetter(3))
        chrm_dict[chrm_n]=sorted(chrm_dict[chrm_n], key=operator.itemgetter(2))
        header=[-1]
        print "sorted"
        for i in range(len(chrm_dict[chrm_n])-1):
            try:
	        dist=chrm_dict[chrm_n][i+1][3]-chrm_dict[chrm_n][i][3]#calc distance between markers
            except:
	        dist=-1
            header.append(dist)
            chrm_dict[chrm_n][i+1].append(dist)
        matri=[]
        for ite in indiv_n: #starts each line with id number
            matri.append([ite])
        snp_order=[]
        for snp in chrm_dict[chrm_n]:#lists snps in order
            snp_order.append(str(snp[0]))
        for ids in matri:
            idss=ids[0]
            for snp in snp_order:
                try:
                    ids.append(snp_dict[(snp,idss)])
                except KeyError: 
                    ids.append('10')
        newdat=matri
        allele1=copy.deepcopy(newdat[:])
        allele2=copy.deepcopy(newdat[:])
        for x in range(len(newdat)):
            for y in range(len(newdat[x])):
                if newdat[x][y]=='1':
                    allele1[x][y]='1'
                    allele2[x][y]='1'
                elif newdat[x][y]=='2':
                    allele1[x][y]='2'
                    allele2[x][y]='2'
                elif newdat[x][y]=='3':
                    allele1[x][y]='1'
                    allele2[x][y]='2'
                elif newdat[x][y]=='10':
                    allele1[x][y]='-9'
                    allele2[x][y]='-9'
                elif newdat[x][y]=='':
                    allele1[x][y]='-9'
                    allele2[x][y]='-9'
        structdat1=[]
        structdat1=allele1+allele2
        structdat=copy.deepcopy(structdat1)
        structdat.sort()
        for x in range(len(structdat)):
            structdat[x].insert(1,structdat[x][0][0:3])
        for x in range(len(structdat)):
            if (int(structdat[x][0])>121488000 and int(structdat[x][0])<121488999):
                structdat[x][1]='999'
        for x in range(len(structdat)):
            if (int(structdat[x][0])>=999999000 and int(structdat[x][0])<999999100):
                structdat[x][1]='998'
        outputfile='%s.txt1'%nm
        if os.path.exists(outputfile): pass
        else:
          structdat.insert(0,copy.deepcopy(snp_order))
          structdat.insert(1, header)        
        fs=open(outputfile, 'a')
        for x in range(len(structdat)):
            for y in range(len(structdat[x])):
                fs.write(str(structdat[x][y]))
                fs.write(' ')
            fs.write('\n')
        fs.close()




#actually makes the structure files

#for infile in glob.glob( os.path.join(path, '*.txt') ):
#    print infile
#    raw_to_struct(infile)


def get_all(infiles):
    snps=set()
    inds=set()
    for fil in infiles:
        for lin in open(fil):
            snp=lin.split(',')[0].strip('[')
            ind=lin.split(',')[1].strip()
            snps.add(snp)
            inds.add(ind)   
    return(snps,inds)



def count_missing(infiles):
    snp_miss={}
    ind_miss={}
    for item in inds: ind_miss[item]=0
    for item in snps: snp_miss[item]=0
    for fil in infiles:
      try:
        print fil
        for lin in open(fil):
            snp=lin.split(',')[0].strip('[')
            ind=lin.split(',')[1].strip()
            geno=int(lin.split(',')[2].strip(']\n'))
            if geno != 10:
	        ind_miss[ind]=ind_miss[ind]+1
		snp_miss[snp]=snp_miss[snp]+1
      except:
        print "EROOR"
        print fil
        print lin
    return (snp_miss,ind_miss)




#inf=[]
#for infile in glob.glob( os.path.join('./RawData_PNAS', '*.txt') ):
#    print infile
#    inf.append(infile)


#count_missing(inf)




#def raw_to_phase_inp(infile,outfile):





#def do_phasing(infile, outfile):




def phased_to_raw(phasedir,chrm_range=range(1,31),k3switch='no'):
    (ibmclocdict, locdict, nam_dict,trans,chrm_dict,chrm_n,forw_dict,rev_dict,inddict,ABdict)=make_dicts()
    revdict3k=make3k_revdict()
    for item in chrm_range:
        oufi=open('%s/Raw_phased_chrmA%i.txt'%(phasedir,item),'w')
	print item
	if k3switch=='yes':
            snps=open("%s/chrm%s_3k.txt1"%(phasedir,item)).readline().split() #opens structure input file inorder to pull SNP locations iCOORDINTAE SYSTEM!!!        
            hg=open("%s/chrm_3k%s_hapguess_switch.out"%(phasedir,item)).readlines() #opens fastphase outputfile  
        else:   
	    snps=open("%s/chrm%s.txt1"%(phasedir,item)).readline().split() #opens structure input file inorder to pull SNP locations iCOORDINTAE SYSTEM!!!        
            hg=open("%s/chrm%s_hapguess_switch.out"%(phasedir,item)).readlines() #opens fastphase outputfile  
        assert len(snps)==len(hg[22].split())#makes sure number of snps in head == number in lines
        for x,lin in enumerate(hg):#skipps all the headers and footers
           if lin.startswith('# id'):
                idi=lin.split()[2]
		for i,snp in enumerate(snps):
                    A=hg[x+1].split()[i]
                    B=hg[x+2].split()[i]
		    geno=ABdict[A,B]
		    if k3switch=='no':
			osnp=rev_dict[snp] #back to icmb numbers from SNPmap
		    if k3switch=='yes':
			osnp=revdict3k[snp]
                    oufi.write("["+osnp+','+idi+','+geno+']'+'\n')
           else:pass             






#phasedir=("/share/HillisLab/ebm447/EIG4.2/POPGEN/55k")






def raw_to_eig(infiles,outstr,group='region'):#infiles is a list of files
    (ibmclocdict, locdict, nam_dict,trans,chrm_dict,chrm_n,forw_dict,rev_dict,inddict,ABdict)=make_dicts()
    if type(infiles) != list:
	print "infiles needs to be alist"
    snps=set()
    ids=set()
    snpfi=open(outstr+'.snp','w')
    goufi=open(outstr+'.geno','w')
    indfi=open(outstr+'.ind','w')
    for infile in infiles:
        goufi=open(outstr+'.geno','a')
	indfi=open(outstr+'.ind','a')
        print infile
        fi=open(infile)
        rawdat=[]
	for lin in fi:
          try:
                stri=str(str(lin)[1:-2])
                stri2=(str(stri)).split(',')
                snpnum=locdict[ibmclocdict[stri2[0].strip()]]#this is trnslating the ibmc SNP number to that from SNPmap
                lst=[snpnum,stri2[1].strip(),trans[stri2[2].strip()]]
                ids.add(stri2[1].strip())
		snps.add(snpnum)
		if (lst[2]!='10'):
                                rawdat.append(lst)
          except:pass
        for lin in rawdat:
                goufi.write(lin[0]+' '+lin[1]+' '+lin[2]+'\n')
        goufi.close()
    for item in ids:
        if group=='region':
		indfi.write(item+'\tF\t'+inddict[item][4].strip('"')+'\n')#by region
        
	if group =='breed':
		indfi.write(item+'\tF\t'+inddict[item][2].strip('"')+'\n')#by breedcode
	if group =='country':
		indfi.write(item+'\tF\t'+inddict[item][5].strip('"')+'\n')
    snp_info=open('SNP_Map.txt','r').readlines()
    for lin in snp_info:
  	lii=lin.split('\t')
  	if lii[0] in snps:
		chrm=lii[2]
		if lii[2]=='"X"': chrm='30'
   		snpfi.write(" ".join([lii[0],chrm,'0.0',lii[3]])+'\n') #pulls outsnp number, chrm, dummy recombination distance, lcoaction
    snpfi.close()
    par=open("par.example",'r')
    opar=open('par.'+outstr,'w')
    for lin in par.readlines():
		opar.write(lin.replace('example',outstr))
    opar.close()


'''#infs=[]
#for infile in glob.glob( os.path.join('./Feb9_3k', 'Raw_ph*.txt') ):
#    print infile
#    infs.append(infile)



#raw_to_eig(infs,'raweig_3k')
'''
'''for item in raw_direc:
    raw_to_phase_inp(item,tmp)
    do_phasing(tmp,out)
'''



def eig_subsamp(instr,outstr,num=30,group='breedcode'):#infiles is a list of files
    (ibmclocdict, locdict, nam_dict,trans,chrm_dict,chrm_n,forw_dict,rev_dict,inddict,ABdict)=make_dicts()
    snpfi=open(instr+'.snp','r')
    goufi=open(instr+'.geno','r')
    indfi=open(instr+'.ind','r')
    par=open('par.'+instr,'r')
    osnpfi=open(outstr+'.snp','w')
    ogoufi=open(outstr+'.geno','w')
    oindfi=open(outstr+'.ind','w')
    opar=open('par.'+outstr,'w')
    if group=='breedcode':ref=2
    if group=='region':ref=4
    if group=='country':ref=5
    groups_dict={}
    subinds=set()
    for lin in indfi:
	ind=lin.split()[0]
	ke=inddict[ind][ref]
        if ke not in groups_dict: groups_dict[ke]=[]
        groups_dict[ke].append(lin)
    for kes in groups_dict:
	print "sampling %s" %kes
	if len(groups_dict[kes])>num:
	    sub=random.sample(groups_dict[kes], num)
	else: sub=groups_dict[kes]
	for item in sub:
	    oindfi.write(item)
	    subinds.add(item.split()[0])
    oindfi.close()
    subsnps=set()
    print "stripping geno file"
    for lii in goufi:
	if lii.split()[1] in subinds:
	    ogoufi.write(lii)
	    subsnps.add(lii.split()[0])
    ogoufi.close()
    print "stripping snp file"
    for lii in snpfi:
	if lii.split()[0]in subsnps:
		osnpfi.write(lii)
    print "writing par file"
    for lin in par.readlines():
                opar.write(lin.replace(instr,outstr))
    opar.close()




def SNP_subsamp(instr,outstr,output,chrm,interval=50,wind_size=50,group='breedcode'):#infiles is a list of files
    (ibmclocdict, locdict, nam_dict,trans,chrm_dict,chrm_n,forw_dict,rev_dict,inddict,ABdict)=make_dicts()
    goufi=open(instr+'.geno','r').readlines()
    indfi=instr+'.ind'
    inds=set()
    indifi=open(indfi).readlines()
    win_dict={}  
    for lin in indifi:
        inds.add(lin.split()[0])
        win_dict[lin.split()[0]]={}
    par=open('par.'+instr,'r').readlines()
    windows=len(chrm_dict[chrm])
    winds=[]
    bps=[]
    for wind in range(0,windows,interval):
      winds.append(wind)
      snpz=chrm_dict[chrm][wind:wind+wind_size]
      bp=snpz[0][3]
      bps.append(bp)
      snps=set()
      osnpfi=open(outstr+'%i.snp'%wind,'w')
      for lin in snpz:
        snps.add(lin[0])
        assert chrm==lin[2]
        if int(chrm) > 22:
          osnpfi.write(" ".join([str(lin[0]),chrm, '0.0',str(lin[3]),'\n']))
        else:
          osnpfi.write(" ".join([str(lin[0]),'1', '0.0',str(lin[3]),'\n']))
      osnpfi.close()
      ogoufi=open(outstr+'%i.geno'%wind,'w')
      for geno in goufi:
         if int(geno.split()[0]) in snps:
           ogoufi.write(geno)    
      ogoufi.close()
      oindfi=outstr+'%i.ind'%wind
      os.system('cp %s %s' %(indfi,oindfi))
      opars='par.'+"".join([outstr+'%i'%wind])
      opar=open(opars,'w')
      for lin in par:
                opar.write(lin.replace(instr,"".join([outstr+'%i'%wind])))            
      opar.close()
      os.system("/home/ejbm/Documents/EIG4.2/bin/smartpca -p %s" %opars)      
      evec=open('%s%i.evec' %(outstr,wind)).readlines()
      indics=[]
      try:
        for lin in evec[1:]:
          indics.append(float(lin.split()[1]))
        ran=max(indics)-min(indics)
        for lin in evec[1:]:
          lii=lin.split()
          win_dict[lii[0]][bp]=(float(lii[1])-min(indics))/ran 
      except: pass
    oufi=open(output,'w')
    for it in bps: oufi.write(', '+str(it))
    oufi.write('\n')   
    for ki in win_dict:
        oufi.write(ki+', ')
        for bp in bps:   
            try:
               oufi.write(str(win_dict[ki][bp])+', ')
            except:
                oufi.write('ERR, ')	      
        oufi.write('\n')
    oufi.close()
    
