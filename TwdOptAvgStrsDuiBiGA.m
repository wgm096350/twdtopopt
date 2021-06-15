clc
clear
addpath ..;
addpath ..\mma;
load twd_final.mat;
global inputstrs % cell for visualization of stress, set it of global scope for conveniency
volscale=1000;strsscale=1200;filterwidth=28;loadedgeupcoord=24;
TaxisOptMass=3.2460E-02;TaxisOptVol=TaxisOptMass/rho;
strsyield=1170;strsfracture=1310;%Yielding Stress
safetyfactor=1.05;poliesafetyfactor=1.67;%Safety Factor
fixloadcols=5;%number of columns fixed at the loading location
% load('twdoptres_212MPa.mat')

stressp=strsyield/safetyfactor;stresspolie=strsfracture/poliesafetyfactor;%maximum allowable Mises stress
penal=1;%penal factor for elastic modulus
penal2=1;%penal factor for material density(ML^-3)
E0=E*1e-3;rho0=rho*1e-3;%the elastic modulus and density for void element
elecolumns=size(mp,2)-fixloadcols;%size(mp,2) is total element columns, 'elecolumns' is the number of active element columns whose width is changeable
elerows=size(designdomain,1)/size(mp,2);%'elerows' is the number of element rows
element_1=reshape(ELEMENT(:,1),size(mp,2),elerows);
% obtain the number of elements whose density is fixed to be 1, ie those
% element locating at loading point, ie non-design domain
nondesigndomain=[];
for jjj=1:fixloadcols
    nondesigndomain=[nondesigndomain,element_1(end-jjj+1,1:length(nset('load1'))-1-0*(jjj-1))];
end
nondesigndomain=nondesigndomain(:);
element_1(end-fixloadcols+1:end,:)=[];designdomain=element_1(:);% extract element Nos in design domain
%Forming filtering-related matrix H and Hs
% fn=31;
% iH=repmat((1:elecolumns)',fn,1);
% jH=bsxfun(@plus,repmat(1:fn,elecolumns,1),(0:elecolumns-1)');
% sH=repmat([1:2:fn,fn-2:-2:1]*0.5,elecolumns,1);
% H=sparse(iH,jH,sH);
% H(:,[1:(fn-1)/2,end-(fn-1)/2+1:end])=[];%delete columns [1:(fn-1)/2,end-(fn-1)/2+1:end] in H
% Hs=sum(H,2);%calculate summation of each row
% fn=21;
% iH=repmat((1:size(mp,2))',fn,1);
% jH=bsxfun(@plus,repmat(1:fn,size(mp,2),1),(0:size(mp,2)-1)');
% sH=repmat([1:2:fn,fn-2:-2:1]*0.5,size(mp,2),1);
% H=sparse(iH,jH,sH);
% H(:,[1:(fn-1)/2,end-(fn-1)/2+1:end])=[];%delete columns [1:(fn-1)/2,end-(fn-1)/2+1:end] in H
% Hs=sum(H,2);%calculate summation of each row
% H=bsxfun(@rdivide,H,Hs);
% nonfixcolumns=1:elecolumns;fixcolumns=elecolumns+1:size(mp,2);
% Hvv=H(nonfixcolumns,nonfixcolumns);Hvf=H(nonfixcolumns,fixcolumns);
% Hfv=H(fixcolumns,nonfixcolumns);Hff=H(fixcolumns,fixcolumns);
% Heff=Hvv-Hvf*inv(Hff)*Hfv;
%helmholtz type filtering
eleL=(max(node(:,1))-min(node(:,1)))/size(mp,2);
Hs=ones(size(mp,2),1);
nonfixcolumns=1:elecolumns;fixcolumns=elecolumns+1:size(mp,2);

Hmatw=helmholtz1d_swd(size(mp,2),eleL,eleL*filterwidth);
Hmatw=bsxfun(@rdivide,Hmatw,Hs);
Hmatwvv=Hmatw(nonfixcolumns,nonfixcolumns);Hmatwvf=Hmatw(nonfixcolumns,fixcolumns);
Hmatwfv=Hmatw(fixcolumns,nonfixcolumns);Hmatwff=Hmatw(fixcolumns,fixcolumns);
Hmatweff=Hmatwvv-Hmatwvf*inv(Hmatwff)*Hmatwfv;

Hvoidw=helmholtz1d_swd(size(mp,2),eleL,eleL*filterwidth);
Hvoidw=bsxfun(@rdivide,Hvoidw,Hs);
Hvoidwvv=Hvoidw(nonfixcolumns,nonfixcolumns);Hvoidwvf=Hvoidw(nonfixcolumns,fixcolumns);
Hvoidwfv=Hvoidw(fixcolumns,nonfixcolumns);Hvoidwff=Hvoidw(fixcolumns,fixcolumns);
Hvoidweff=Hvoidwvv-Hvoidwvf*inv(Hvoidwff)*Hvoidwfv;
%construct cell 'input' for visualization
input.node=node;
input.element=element;
input.Output_DensityField=1;
input.isovalue=0;
input.Output_DeformationField=0;
inputstrs=input;
%the first variablenum values in xmax and xmin denotes the upper/lower
%bound of each column's material width;the second variablenum values in xmax
%and xmin denotes the upper/lower bound of each column's inner void width
rampwidth=1/elerows;%width of the 2 ramps,typically equals to 1/(element rows)
xmax=[0.68*ones(elecolumns,1);0.3*ones(elecolumns,1)];
xmin=[1/elerows*ones(elecolumns,1);0.0*ones(elecolumns,1)];
low=xmin;upp=xmax;%asymptote limit for MMA procudure
%initialization of design variale,'mat_w' denotes the material width of each
%column while 'invoid_w' denotes the width of each column's inner void

mat_w=1/elerows+(0.60-1/elerows)*rand(elecolumns,1);mat_initial=mat_w;
invoid_w=0.00+(0.30-0.00)*rand(elecolumns,1);invoid_initial=invoid_w;
load 'twdoptres_212MPa.mat' mat_initial invoid_initial
mat_w=mat_initial;invoid_w=invoid_initial;
% load senstest.mat
% mat_w(45)=mat_w(45)+0.5e-6;
% invoid_w(145)=invoid_w(145)+0.5e-6;
%obtain the material height of each column (excluding fixed columns) by interpolating the height from
%optimization result of T-Axis. Assumming each column's width is equal to eleL.
% middlecolumnR=(min(node(:,1))+eleL/2):eleL:(max(node(:,1))-fixloadcols*eleL-eleL/2);
% mat_wfilter_Taxisopt=nihe_cailiaogaodu(middlecolumnR',max(node(:,2)));
% mat_woriginal_Taxisopt=H\(mat_wfilter_Taxisopt.*Hs);
% mat_w=mat_woriginal_Taxisopt;
%filter of design variables to achieve smooth transition
mat_wfixed=loadedgeupcoord/elerows*ones(fixloadcols,1);
invoid_wfixed=0*ones(fixloadcols,1);
% ttttt=inv(Hff)*(mat_wfixed-Hfv*mat_w); % ttttt=5.469
mat_wfilter=Hmatweff*mat_w+Hmatwvf*inv(Hmatwff)*mat_wfixed;
% mat_wfilter=mat_w;
invoid_wfilter=Hvoidweff*invoid_w+Hvoidwvf*inv(Hvoidwff)*invoid_wfixed;

%normalized coordinated of each element, the i-th row of 'normcoord'
%denotes the normalized coordinated of the i-th column of elements,the
%first column is the column closest to the rotational axis
normcoord=reshape(mp*ones(size(mp,2),1),size(mp,2),elerows);
normcoord(end-fixloadcols+1:end,:)=[];
%Extend vector 'mat_wfilter' ,'invoid_wfilter' to the same size of
%normcoord. This extending operation facilitates the vectorization of function 'DensDistributionFunc'
extmat_wfilter=repmat(mat_wfilter,1,elerows);
extinvoid_wfilter=repmat(invoid_wfilter,1,elerows);
%vector of element pseudo density 'xPhys', initialized to 1.0
xPhys=zeros(size(element,1),1);xPhys(nondesigndomain)=1.0;
%Obtain each element's pseudo density saved as 'xx2D' in matrix form
xx2D=DensDistributionFuncFillet(extinvoid_wfilter,rampwidth,normcoord,extmat_wfilter);
% xx2D=DensDistributionFunc(extinvoid_wfilter,rampwidth,normcoord,extmat_wfilter);
%xx2D saves elements' pseudo density in matrix form, the i-th row of it
%denotes the pseudo density of the i-th column of elements. But the element
%numbering sequence is by row. So rearrangement of xx2D is needed.
xPhys(designdomain) = xx2D(:);
% load('swdoptres.mat');
% E0=E*1e-6;rho0=rho*1e-6;
% xPhys(:) = 1;
%visualization of element pseudo density distribution
input.x=xPhys;
figure(1);
postprocessing(input);

%construct vector of design variables 'xold1','xold2' equals to 'xold1'
%initially, both used in MMA solver;'change' is the convergence tolerance
xold1=[mat_w;invoid_w];
xold2=xold1;
change=1;

loop=0;%iteration indicator
Compliance=[];
rholist=[];
xlist=[];
filteredmatvoidw=[];
nofilteredmatvoidw=[];
while change > 0.001 && loop<200
    loop = loop+1;
    filteredmatvoidw=[filteredmatvoidw,[mat_wfilter;invoid_wfilter]];
    nofilteredmatvoidw=[nofilteredmatvoidw,[mat_w;invoid_w]];
    xlist=[xlist,xPhys(:)];
    
%     load('twdoptres_212MPa.mat')
%     xPhys=xlist(:,end);
%     critrho=0.5;
%     xPhys(xPhys>critrho)=1;
%     xPhys(xPhys<critrho)=0;

    [iK,jK,sK]=assemble_AXQ4_mex(element,node,E0+(E-E0)*xPhys.^penal,nu*ones(size(xPhys)));
    K=sparse(iK,jK,sK);K=(K+K')/2;
    [iF,jF,sF,V]=assemble_AXQ4_inertia_mex(element,node,rho0+(rho-rho0)*xPhys.^penal2,w);
    force_inertia=sparse(iF,jF,sF,size(node,1)*2,1);
    F=force1+force_inertia;
    U=zeros(size(F));U(freedofs)=K(freedofs,freedofs)\F(freedofs,1);
    %     Compliance(loop)=full(U'*K*U)/objscale;
    VolumeCons(loop)=xPhys'*V/TaxisOptVol*volscale;
    
    figure(2);
    subplot(2,1,1)
    plot(1:loop,VolumeCons)
    title('VolumeCons');
    
    dkr=penal*(E-E0)*xPhys.^(penal-1)./(E0+(E-E0)*xPhys.^penal);
    dfr=penal2*(rho-rho0)*xPhys.^(penal2-1)./(rho0+(rho-rho0)*xPhys.^penal2);
    
    %     ComplianceSensitivity=sensitivity_AXQ4_inertia(U,edofMat,sF,dfr)+sensitivity_AXQ4(U,edofMat,sK,dkr);
    %     ComplianceSensitivity=ComplianceSensitivity(designdomain)/objscale;
    VolumeSensitivty=V(designdomain)/TaxisOptVol*volscale;
    
    %     liebiao=304:elecolumns:(304+elecolumns*(elerows-1));
    %     sensilist_lastcolumn(:,loop)=ComplianceSensitivity(liebiao);
    %     xPhys_lastcolumn(:,loop)=xPhys(liebiao);
    
    q=0;p=8;
    dE=E*(penal-q)*xPhys.^(penal-q-1);
    %compute normalized p-norm Mises stress 'maxstress', real maximum Mises
    %stress 'maxstress_exact'(also normalized), 'dmaxstress' is sensitivity of 'maxstress'
    [pnormMaxStrs,dpnormMaxStrsdxPhys,maxstrs_exact,avghoop,davghoopdxPhys,pnormMaxAvgSigrr,dpnormMaxAvgSigrrPhys,maxavgsigrr_exact]=maxstress_calc_AXQ4_4_int(xPhys,E*xPhys.^(penal-q),...
        dE,nu*ones(size(xPhys)),node,element,U,K,sK,sF,edofMat,freedofs,dkr,dfr,p,stressp,stresspolie,1:size(element,1),size(mp,2),size(element,1)/size(mp,2));
    k=maxstrs_exact/pnormMaxStrs;
%     k=1;
    StrsCons(loop)=(pnormMaxStrs*k-1)*strsscale;%This actually equals to 'cs=(maxstress_exact-1)*100;'
    
    AvgHoopStrsCons(loop)=(avghoop-1)*strsscale;%This actually equals to 'cs=(maxstress_exact-1)*100;'
    
    krr=maxavgsigrr_exact/pnormMaxAvgSigrr;
%     krr=1;
    MaxAvgSigrrCons(loop)=(pnormMaxAvgSigrr*krr-1)*strsscale;
    
    figure(2);
    subplot(2,1,2)
    plot(1:loop,StrsCons)
    title('StrsCons');
    
    %Since the sensitivity of 'maxstress_exact' can't be caluculated, so
    %it's approximated. 'dmaxstress' is the approximated sensitivity of 'maxstress_exact'
    dmaxstrsdxPhys=dpnormMaxStrsdxPhys(designdomain)*k*strsscale;
    davghoopstrsdxPhys=davghoopdxPhys(designdomain)*strsscale;
    dmaxavgsigrrdxPhys=dpnormMaxAvgSigrrPhys(designdomain)*krr*strsscale;
    %construct vector of constraints
    ConsVec=[StrsCons(loop) AvgHoopStrsCons(loop) MaxAvgSigrrCons(loop)]';
    %Construciton and vectorization of the matrix consisting of element's
    %sensitivity w.r.t inner void width (shiftdiff2) and material width (shiftdiff)
    shiftdiff_iF=1:1:size(designdomain,1);
    shiftdiff_jF=repmat(1:1:elecolumns,1,elerows);
    dxPhysdmat_w=DensDistributionFuncFilletDiffMatWidth(extinvoid_wfilter,rampwidth,normcoord,extmat_wfilter);dxPhysdmat_w=dxPhysdmat_w(:);
    %     dxPhysdmat_w=DensDistributionFuncDiffMatWidth(extinvoid_wfilter,rampwidth,normcoord,extmat_wfilter);dxPhysdmat_w=dxPhysdmat_w(:);
    dxPhysdmat_w=sparse(shiftdiff_iF,shiftdiff_jF,dxPhysdmat_w);
    dxPhysdinvoid_w=DensDistributionFuncFilletDiffVoidWidth(extinvoid_wfilter,rampwidth,normcoord,extmat_wfilter);dxPhysdinvoid_w=dxPhysdinvoid_w(:);
    %     dxPhysdmat_w=DensDistributionFuncDiffMatWidth(extinvoid_wfilter,rampwidth,normcoord,extmat_wfilter);dxPhysdmat_w=dxPhysdmat_w(:);
    dxPhysdinvoid_w=sparse(shiftdiff_iF,shiftdiff_jF,dxPhysdinvoid_w);
    
    %transform sensitivity w.r.t xPhys to design variables
    %     dcompliancedmat_w=(ComplianceSensitivity'*dxPhysdmat_w)';
    %     dcompliancedmat_w(:)=Heff*(dcompliancedmat_w(:));
    %     ComplianceSensitivity=dcompliancedmat_w;
    
    %     sensilist_5(:,loop)=ComplianceSensitivity(end-4:end);
    
    dmaxstrsdmat_w=(dmaxstrsdxPhys'*dxPhysdmat_w)';
    dmaxstrsdmat_w(:)=Hmatweff*(dmaxstrsdmat_w(:));
    dmaxstrsdinvoid_w=(dmaxstrsdxPhys'*dxPhysdinvoid_w)';
    dmaxstrsdinvoid_w(:)=Hvoidweff*(dmaxstrsdinvoid_w(:));
    MaxstrsSensitivity=[dmaxstrsdmat_w;dmaxstrsdinvoid_w];
    
    davgtrsdmat_w=(davghoopstrsdxPhys'*dxPhysdmat_w)';
    davgtrsdmat_w(:)=Hmatweff*(davgtrsdmat_w(:));
    davgstrsdinvoid_w=(davghoopstrsdxPhys'*dxPhysdinvoid_w)';
    davgstrsdinvoid_w(:)=Hvoidweff*(davgstrsdinvoid_w(:));
    AvgstrsSensitivity=[davgtrsdmat_w;davgstrsdinvoid_w];
    %     MaxstrsSensitivity=[dmaxstrsdmat_w;0*ones(elecolumns,1)];
    
    dmaxavgsigrrdmat_w=(dmaxavgsigrrdxPhys'*dxPhysdmat_w)';
    dmaxavgsigrrdmat_w(:)=Hmatweff*(dmaxavgsigrrdmat_w(:));
    dmaxavgsigrrdinvoid_w=(dmaxavgsigrrdxPhys'*dxPhysdinvoid_w)';
    dmaxavgsigrrdinvoid_w(:)=Hvoidweff*(dmaxavgsigrrdinvoid_w(:));
    MaxAvgSigrrSensitivity=[dmaxavgsigrrdmat_w;dmaxavgsigrrdinvoid_w];
    
    dVdmat_w=(VolumeSensitivty(:)'*dxPhysdmat_w)';
    dVdmat_w(:)=Hmatweff*(dVdmat_w(:));
    dVdinvoid_w=(VolumeSensitivty(:)'*dxPhysdinvoid_w)';
    dVdinvoid_w(:)=Hvoidweff*(dVdinvoid_w(:));
    VolumeSensitivty=[dVdmat_w;dVdinvoid_w];
    %     VolumeSensitivty=[dVdmat_w;0*ones(elecolumns,1)];
    
    
    n=2*elecolumns;%n denotes number of design variables
    dConsVec=[MaxstrsSensitivity(:)';AvgstrsSensitivity';MaxAvgSigrrSensitivity'];
    m=3;%m denotes number of constraint
    a0 = 1;a = zeros(m,1);ccc=1000*ones(m,1);d = ones(m,1);%constant parameters for MMA solver
    
    [xnew,ymma,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(m,n,loop,[mat_w;invoid_w],xmin,xmax,xold1,xold2, ...
        VolumeCons(loop),VolumeSensitivty,ConsVec,dConsVec,low,upp,a0,a,ccc,d);
    xold2=xold1;xold1=[mat_w;invoid_w];
    %update design variables
    mat_w=xnew(1:elecolumns);invoid_w=xnew(elecolumns+1:end);
    %start next iteration
    mat_wfilter=Hmatweff*mat_w+Hmatwvf*inv(Hmatwff)*mat_wfixed;
    invoid_wfilter=Hvoidweff*invoid_w+Hvoidwvf*inv(Hvoidwff)*invoid_wfixed;
    extmat_wfilter=repmat(mat_wfilter,1,elerows);
    extinvoid_wfilter=repmat(invoid_wfilter,1,elerows);
    xx2D=DensDistributionFuncFillet(extinvoid_wfilter,rampwidth,normcoord,extmat_wfilter);
    %     xx2D=DensDistributionFunc(extinvoid_wfilter,rampwidth,normcoord,extmat_wfilter);
    %visualize new density distribution
    xPhys(designdomain) = xx2D(:);
    input.x=xPhys;
    figure(1);
    postprocessing(input);
    hold on;
    plot(node(nset('load1'),1),node(nset('load1'),2),'ro','MarkerSize',2);hold off;
%     figure(5)
%     plot(MaxstrsSensitivity(1:elecolumns))
%     figure(6)
%     plot(MaxstrsSensitivity(elecolumns+1:end))
    %display the location of loaded node
    disp([' It.: ' sprintf('%4i',loop) ' VolCons: ' sprintf('%10.4f',VolumeCons(loop)/volscale)...
        ' StrsCons: ' sprintf('%10.4f',StrsCons(loop)/strsscale)...
        ' HoopStrsCons: ' sprintf('%10.4f',AvgHoopStrsCons(loop)/strsscale)...
        ' RadStrsCons: ' sprintf('%10.4f',MaxAvgSigrrCons(loop)/strsscale)...
        ' ch.: ' sprintf('%6.3f',change )]);
    %update convergence tolerance
    change=max(abs(xnew-xold1));
end
save twdoptres.mat;