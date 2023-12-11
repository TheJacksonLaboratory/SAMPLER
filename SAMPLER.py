import numpy as np
import os
import pickle
from scipy.stats import norm


###################################################################
###helper functions
#####################################################
#####################################################
def SAMPLERrep(X,pers=[5,15,25,35,45,55,65,75,85,95]):
    wsi_SAMPLER_rep=np.ravel(np.percentile(X,pers,axis=0))
    
    return wsi_SAMPLER_rep


def Fscore(feats,pers,Ftab,Stab,Wvec):
    
    pers=np.asarray(pers)
    pers=pers/100
    
    pernum=len(pers)
    tnum=len(feats)
    places=np.argsort(feats)
    ranks=np.arange(tnum)
    
    wsi_pers=np.zeros((tnum,))
    wsi_pers[places]=ranks
    wsi_pers=(wsi_pers+1)/(tnum+2)
    
    scores=np.zeros((tnum,))
    for i in range(tnum):
        
        cfeat=feats[i]
        cper=wsi_pers[i]
        
        
        if cper<np.min(pers)+0.0001:
            scores[i]=Wvec[0]*np.interp(cfeat,Ftab[0,:],Stab[0,:])
        elif cper>np.max(pers)-0.0001:
            scores[i]=Wvec[pernum-1]*np.interp(cfeat,Ftab[pernum-1,:],Stab[pernum-1,:])
        else:
            
            ##find percentile, basic code, issues with indexing and ...
            check=1
            for i in range(pernum):
                if (cper<pers[i]) and (check>0):
                    iplace=i-1
                    check=0
            
            
            iper=pers[iplace]
            ip1per=pers[iplace+1]
            landa=(cper-iper)/(ip1per-iper)


            scoreb=np.interp(cfeat,Ftab[iplace,:],Stab[iplace,:])
            scorea=np.interp(cfeat,Ftab[iplace+1,:],Stab[iplace+1,:])
            
            cL=landa*scoreb+(1-landa)*scorea
            cW=landa*Wvec[iplace]+(1-landa)*Wvec[iplace+1]
            
            scores[i]=cW*cL
    
    return scores


def Sscore(Xmat,Xtab,Ltab,Wmat,pers):
    
    tnum, fnum =np.shape(Xmat)
    
    scores=np.zeros((tnum,fnum))
    
    for i in range(fnum):
        feats=Xmat[:,i]
        Ftab=np.squeeze(Xtab[i,:,:])
        Stab=np.squeeze(Ltab[i,:,:])
        Wvec=Wmat[:,i]
        scores[:,i]=Fscore(feats,pers,Ftab,Stab,Wvec)
    
    return scores


def data_loader(spath,ddim):
    
    tiles=os.listdir(spath)
    tnum=len(tiles)
    
    XR=np.zeros((tnum,ddim),dtype=np.float32)
    XP=np.zeros((tnum,ddim),dtype=np.float32)
    XB=np.zeros((tnum,ddim),dtype=np.float32)
    
    locs=np.zeros((tnum,2),dtype=np.int32)
    
    for i in range(tnum):
        
        tile=tiles[i]
        
        W=pickle.load(open(spath+tile,'rb'))
        XR[i,:]=W['feats'][0,:]
        XP[i,:]=W['feats'][1,:]
        XB[i,:]=W['feats'][2,:]
        
        locs[i,:]=W['locs']
    
    return XR, XP, XB, locs


def hmapgen3(scores,locs,Tsize,center=True):
    
    if 1*center>0:
        bias=Tsize
    else:
        bias=0
    
    nr,nc=np.max(locs,axis=0)+Tsize+2*bias
    
    hmap3=np.zeros((nr,nc,3))
    nmap=np.zeros((nr,nc))
    
    for i in range(len(locs)):
        ci=locs[i,0]+bias
        cj=locs[i,1]+bias
        
        
        hmap3[ci:ci+Tsize,cj:cj+Tsize,0]=hmap3[ci:ci+Tsize,cj:cj+Tsize,0]+scores[i]
        hmap3[ci:ci+Tsize,cj:cj+Tsize,2]=hmap3[ci:ci+Tsize,cj:cj+Tsize,2]+1-scores[i]
        nmap[ci:ci+Tsize,cj:cj+Tsize]=nmap[ci:ci+Tsize,cj:cj+Tsize]+1
    
    hmap3=hmap3/(nmap[:,:,np.newaxis]+0.000000000001)
    for i in range(3):
        cmap=hmap3[:,:,i]
        cmap[nmap<0.5]=1
        hmap3[:,:,i]=cmap
    
    hmap3[hmap3>1]=1
    hmap3[hmap3<0]=0
    
    return hmap3




def normalize_scores(scores,pars=[5,5]):
    
    bottom=pars[0]
    top=100-pars[1]
    scores=scores-np.percentile(scores,bottom)
    scores[scores<0]=0
    scores=scores/np.percentile(scores,top)
    scores[scores>1]=1

    
    return scores



def LLRgen2class_singlescale(Xtrain,Ytrain,fdim,pers,Nres=100,bws=0.1):
    
    
    pernum=len(pers)
    REPS=1  ###regularization value
    
    Xtab=np.zeros((fdim,pernum,Nres),dtype=np.float32)
    Ltab=np.zeros((fdim,pernum,Nres),dtype=np.float32)


    for i in range(fdim):
        for p in range(pernum):
        
            cpone=p*fdim+i
    
            cXtrain=Xtrain[:,cpone]
        
            cstdR=np.std(cXtrain)
    
            cx0=cXtrain[Ytrain==0]
            cx1=cXtrain[Ytrain==1]
        
            ###############################
            ###table of red
            Xtab[i,p,:]=np.linspace(0,np.max(cXtrain)+0.1,Nres)
            for j in range(Nres):
                cx=Xtab[i,p,j]
                L0=np.mean(norm.pdf(cx0,loc=cx,scale=bws*cstdR))+REPS
                L1=np.mean(norm.pdf(cx1,loc=cx,scale=bws*cstdR))+REPS
                Ltab[i,p,j]=np.log(L1/L0)
    
    Dtab={}
    Dtab['Xtab']=Xtab
    Dtab['Ltab']=Ltab
    
    return Dtab



def LLRgen_singlescale(Xtrain,Ytrain,fdim,pers,Nres=100,bws=0.1):
    
    pernum=len(pers)
    numclasses=int(np.max(Ytrain)+1)
    REPS=1  ###regularization value
    
    Dtab_base={}
    
    
    ######################
    ##compute likelihoods
    for cclass in range(numclasses):
        
    
        Xtab=np.zeros((fdim,pernum,Nres),dtype=np.float32)
        Ltab=np.zeros((fdim,pernum,Nres),dtype=np.float32)
        tot_Ltab=np.zeros((fdim,pernum,Nres),dtype=np.float32)
        
        for i in range(fdim):
            for p in range(pernum):
            
                cpone=p*fdim+i
        
                cXtrain=Xtrain[:,cpone]
            
                cstdR=np.std(cXtrain)
        
                cx=cXtrain[Ytrain==cclass,:]
            
                ###############################
                ###table of red
                Xtab[i,p,:]=np.linspace(0,np.max(cXtrain)+0.1,Nres)
                for j in range(Nres):
                    cx=Xtab[i,p,j]
                    LC=np.mean(norm.pdf(cx,loc=cx,scale=bws*cstdR))+REPS
                    Ltab[i,p,j]=LC
        
        tot_Ltab=tot_Ltab+Ltab
        Dtab_base['Xtab']=Xtab
        Dtab_base['LC_'+str(cclass)]=Ltab
        Dtab_base['tot_Ltab']=tot_Ltab
    
    
    #########################
    ##now convert to likelihood ratios
    Dtab={}
    Dtab['Xtab']=Xtab
    tot_Ltab=Dtab_base['tot_Ltab']
    for cclass in range(numclasses):
        Lother_classes=(tot_Ltab-Dtab_base['LC_'+str(cclass)])/(numclasses-1)
        Dtab['LLR_'+str(cclass)]=np.log(Dtab_base['LC_'+str(cclass)]/Lother_classes)
    
    
    return Dtab


def LLRgen_multiscale(Xtrain,Ytrain,P,fdim,pers,Nres=100,bws=0.1,PT=0.001):
    
    pernum=len(pers)
    REPS=1  ###regularization value
    sdim=fdim*pernum
    numclasses=(np.max(Ytrain)+1)
    
    Dtab_base={}
    
    for cclass in range(numclasses):
        ##now generate tile level scores look up table
        XtabR=np.zeros((fdim,pernum,Nres),dtype=np.float32)
        LtabR=np.ones((fdim,pernum,Nres),dtype=np.float32)
        tot_LtabR=np.zeros((fdim,pernum,Nres),dtype=np.float32)
        
        XtabP=np.zeros((fdim,pernum,Nres),dtype=np.float32)
        LtabP=np.ones((fdim,pernum,Nres),dtype=np.float32)
        tot_LtabP=np.zeros((fdim,pernum,Nres),dtype=np.float32)
        
        XtabB=np.zeros((fdim,pernum,Nres),dtype=np.float32)
        LtabB=np.ones((fdim,pernum,Nres),dtype=np.float32)
        tot_LtabB=np.zeros((fdim,pernum,Nres),dtype=np.float32)
        
        for i in range(fdim):
            for p in range(pernum):
                
                cpone=p*fdim+i
            
                cXtrainR=Xtrain[:,cpone]
                cXtrainP=Xtrain[:,sdim+cpone]
                cXtrainB=Xtrain[:,2*sdim+cpone]
                
                cstdR=np.std(cXtrainR)
                cstdP=np.std(cXtrainP)
                cstdB=np.std(cXtrainB)
            
                cxR=cXtrainR[Ytrain==cclass]
                
                cxP=cXtrainP[Ytrain==cclass]
                
                cxB=cXtrainB[Ytrain==cclass]
                
                ###############################
                ###table of red
                XtabR[i,p,:]=np.linspace(0,np.max(cXtrainR)+0.1,Nres)
                if P[i]<PT:
                    for j in range(Nres):
                        cx=XtabR[i,p,j]
                        LC=np.mean(norm.pdf(cxR,loc=cx,scale=bws*cstdR))+REPS
                        LtabR[i,p,j]=LC
                    
                ###############################
                ###table of purple
                XtabP[i,p,:]=np.linspace(0,np.max(cXtrainP)+0.1,Nres)
                if P[sdim+i]<PT:
                    for j in range(Nres):
                        cx=XtabP[i,p,j]
                        LC=np.mean(norm.pdf(cxP,loc=cx,scale=bws*cstdP))+REPS
                        LtabP[i,p,j]=LC
                
                ###############################
                ###table of black
                XtabB[i,p,:]=np.linspace(0,np.max(cXtrainB)+0.1,Nres)
                if P[2*sdim+i]<PT:
                    for j in range(Nres):
                        cx=XtabB[i,p,j]
                        LC=np.mean(norm.pdf(cxB,loc=cx,scale=bws*cstdB))+REPS
                        LtabB[i,p,j]=LC
                    
        tot_LtabR=tot_LtabR+LtabR
        tot_LtabP=tot_LtabP+LtabP
        tot_LtabB=tot_LtabB+LtabB
            
        Dtab_base['XtabR']=XtabR
        Dtab_base['LtabRC_'+str(cclass)]=LtabR
        Dtab_base['tot_LtabR']=tot_LtabR
        
        Dtab_base['XtabP']=XtabP
        Dtab_base['LtabPC_'+str(cclass)]=LtabP
        Dtab_base['tot_LtabP']=tot_LtabP
        
        Dtab_base['XtabB']=XtabB
        Dtab_base['LtabBC_'+str(cclass)]=LtabB
        Dtab_base['tot_LtabB']=tot_LtabB      
        
    #################################
    ###now clean
    Dtab={}
    Dtab['XtabR']=Dtab_base['XtabR']
    Dtab['XtabP']=Dtab_base['XtabP']
    Dtab['XtabB']=Dtab_base['XtabB']
    
    tot_LtabR=Dtab_base['tot_LtabR']
    tot_LtabP=Dtab_base['tot_LtabP']
    tot_LtabB=Dtab_base['tot_LtabB']
    
    
    for cclass in range(numclasses):
        Lother_classesR=(tot_LtabR-Dtab_base['LtabRC_'+str(cclass)])/(numclasses-1)
        LLR=np.log(Dtab_base['LtabRC_'+str(cclass)]/Lother_classesR)
        LLR[np.isnan(LLR)]=0
        LLR[LLR>100]=100
        LLR[LLR<-100]=-100
        Dtab['LLR_R_C_'+str(cclass)]=LLR
        
        
        Lother_classesP=(tot_LtabP-Dtab_base['LtabPC_'+str(cclass)])/(numclasses-1)
        LLR=np.log(Dtab_base['LtabPC_'+str(cclass)]/Lother_classesP)
        LLR[np.isnan(LLR)]=0
        LLR[LLR>100]=100
        LLR[LLR<-100]=-100
        Dtab['LLR_P_C_'+str(cclass)]=LLR
        
        Lother_classesB=(tot_LtabB-Dtab_base['LtabBC_'+str(cclass)])/(numclasses-1)
        LLR=np.log(Dtab_base['LtabBC_'+str(cclass)]/Lother_classesB)
        LLR[np.isnan(LLR)]=0
        LLR[LLR>100]=100
        LLR[LLR<-100]=-100
        Dtab['LLR_B_C_'+str(cclass)]=LLR
    
    return Dtab
