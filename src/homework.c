#include "fem.h"




void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    
    int ierr;
    double r = w/4;
    int idRect = gmshModelOccAddRectangle(0.0,0.0,0.0,w,h,-1,0.0,&ierr); 
    int idDisk = gmshModelOccAddDisk(w/2.0,h/2.0,0.0,r,r,-1,NULL,0,NULL,0,&ierr); 
    int idSlit = gmshModelOccAddRectangle(w/2.0,h/2.0-r,0.0,w,2.0*r,-1,0.0,&ierr); 
    int rect[] = {2,idRect};
    int disk[] = {2,idDisk};
    int slit[] = {2,idSlit};

    gmshModelOccCut(rect,2,disk,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccCut(rect,2,slit,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccSynchronize(&ierr); 

    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
 
    return;
}


double *femElasticitySolve(femProblem *theProblem)
{

    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    
    
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    
    int nLocal = theMesh->nLocalNode;

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    
  //
  //  A faire :-)
  //
    
    //Construction du système linéaire
    for (iElem = 0; iElem < theMesh->nElem; iElem++){
        for (int j = 0; j < nLocal; j++){
            map[j] = theMesh->elem[nLocal*iElem+j];
            x[j] = theNodes->X[map[j]];
            y[j] = theNodes->Y[map[j]];
        }
        
        for (iInteg = 0; iInteg < theRule->n; iInteg++){
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            
            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);
            
            //Interpolation par les fonctions de forme
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0;
            double dydeta = 0;
            
            for (int k = 0; k < theSpace->n; k++){
                dxdxsi += x[k]*dphidxsi[k];
                dxdeta += x[k]*dphideta[k];
                dydxsi += y[k]*dphidxsi[k];
                dydeta += y[k]*dphideta[k];
            }
            
            double jacobian = fabs(dxdxsi*dydeta - dxdeta*dydxsi);
            if (jacobian < 0){   //si jacobian < 0 ça veut dire que nos noeuds sont mal orientés => on réoriente
                int node = theMesh->elem[nLocal*iElem];
                theMesh->elem[nLocal*iElem] = theMesh->elem[nLocal*iElem+2];
                theMesh->elem[nLocal*iElem+2] = node;
            }
            
            for (int i = 0; i < theSpace->n; i++){
                dphidx[i] = (dphidxsi[i]*dydeta - dphideta[i]*dydxsi) / jacobian;
                dphidy[i] = (dphideta[i]*dxdxsi - dphidxsi[i]*dxdeta) / jacobian;
            }
            
            //Assemblage des matrices A et B
            for (int i = 0; i < theSpace->n; i++){
                for (int j = 0; j < theSpace->n; j++){
                    A[2*map[i]][2*map[j]] += weight * jacobian * (a*dphidx[i]*dphidx[j] + c*dphidy[i]*dphidy[j]);
                    A[2*map[i]][2*map[j]+1] += weight * jacobian * (b*dphidx[i]*dphidy[j] + c*dphidy[i]*dphidx[j]);
                    A[2*map[i]+1][2*map[j]] += weight * jacobian * (b*dphidy[i]*dphidx[j] + c*dphidx[i]*dphidy[j]);
                    A[2*map[i]+1][2*map[j]+1] += weight * jacobian * (a*dphidy[i]*dphidy[j] + c*dphidx[i]*dphidx[j]);
                }
                B[2*map[i] + 1] -= weight * jacobian * phi[i] * rho * g; //Pas de force le long de x donc une composante sur deux est nulle ET force de gravité est dirigée vers les y négatifs
            }
            
            for (i = 0; i < theSpace->n; i++) {
                            theSystem->B[2*map[i] + 1] -= weight * jacobian * phi[i] * g * rho;
                            for (j = 0; j < theSpace->n; j++) {
                                theSystem->A[2*map[i]]    [2*map[j]]     += weight*jacobian * (dphidx[i]*a*dphidx[j] + dphidy[i]*c*dphidy[j]);
                                theSystem->A[2*map[i] + 1][2*map[j]]     += weight*jacobian * (dphidy[i]*b*dphidx[j] + dphidx[i]*c*dphidy[j]);
                                theSystem->A[2*map[i]]    [2*map[j] + 1] += weight*jacobian * (dphidx[i]*b*dphidy[j] + dphidy[i]*c*dphidx[j]);
                                theSystem->A[2*map[i] + 1][2*map[j] + 1] += weight*jacobian * (dphidy[i]*a*dphidy[j] + dphidx[i]*c*dphidx[j]);
                            }
            }
        }
                
    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value); }}
    }
    return femFullSystemEliminate(theSystem);

}
