#include "fem.h"

//On peut utiliser gmsh pour le maillage mais pas la résolution de systèmes linéaires

//Notre problème:
//- Un écoduc sur lequel passe un éléphant pour visualiser les déformations
//- On teste plusieurs matériaux avec des modules de young différents
//- On sélectionne le matériau le plus résistant (le moins déformé) et le plus écologique (car ÉCOduc)

//À rajouter:
//Conditions essentielles (Dirichlet) en normale/tengentielle: 1 seul déplacement normal ou tangentiel
//Déplacement normal = 0: équation correspondante : nx * ux + ny *uy = 0 -> au lieu d'avoir 1 dans la matrice A on aura nx et ny (càd sin et cos)

//Combili pour que dans la matrice on ait: nx * ligne x + ny * ligne y
//Tangentielle:
//Quand on impose le déplacement sur une des 2 équations: ex dépl normal = 0

//nx * équation en x + ny * équation en y -> conditions en normale
//-ny * équation en x + nx * équation en y (vérifier les signes) -> conditions en tangentielle
//nx et ny varient pour chaque noeuds = moyenne des normales des deux segments adjacents, pondérée par la longueur des segments
//Exemple: poutre dans le mur: endroit dans le mur = dirichlet à 0, si un poids sur la poutre = Neumann non homogène, partt ailleurs = Neumann homogène
//Conditions normales (Neumann) en xy et en normale/tangentielle

//1) Remplacer nouvelles lignes:  0 0 0 cos,sin 0 0 0 = déplacement normal (tangentiel = -sin cos)
//2) Imposer déplacement normal ou tangentiel en fonction de ce qu'on veut imposer

//modifier principalement les boundary conditions (pas le code en lui-même)



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
   
   
   for (iElem = 0; iElem < theMesh->nElem; iElem++) {
       for (j=0; j < nLocal; j++) {
           map[j]  = theMesh->elem[iElem*nLocal+j];
           mapX[j] = 2*map[j];
           mapY[j] = 2*map[j] + 1;
           x[j]    = theNodes->X[map[j]];
           y[j]    = theNodes->Y[map[j]];}
       
       for (iInteg=0; iInteg < theRule->n; iInteg++) {
           double xsi    = theRule->xsi[iInteg];
           double eta    = theRule->eta[iInteg];
           double weight = theRule->weight[iInteg];
           femDiscretePhi2(theSpace,xsi,eta,phi);
           femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
           
           double dxdxsi = 0.0;
           double dxdeta = 0.0;
           double dydxsi = 0.0;
           double dydeta = 0.0;
           for (i = 0; i < theSpace->n; i++) {
               dxdxsi += x[i]*dphidxsi[i];
               dxdeta += x[i]*dphideta[i];
               dydxsi += y[i]*dphidxsi[i];
               dydeta += y[i]*dphideta[i]; }
           double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
           
           for (i = 0; i < theSpace->n; i++) {
               dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
               dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }
           for (i = 0; i < theSpace->n; i++) {
               for(j = 0; j < theSpace->n; j++) {
                   A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] +
                                           dphidy[i] * c * dphidy[j]) * jac * weight;
                   A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] +
                                           dphidy[i] * c * dphidx[j]) * jac * weight;
                   A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] +
                                           dphidx[i] * c * dphidy[j]) * jac * weight;
                   A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] +
                                           dphidx[i] * c * dphidx[j]) * jac * weight; }}
            for (i = 0; i < theSpace->n; i++) {
               B[mapY[i]] -= phi[i] * g * rho * jac * weight; }}}
 
   int *theConstrainedNodes = theProblem->constrainedNodes;
   for (int i=0; i < theSystem->size; i++) {
       if (theConstrainedNodes[i] != -1) {
           double value = theProblem->conditions[theConstrainedNodes[i]]->value;
           femFullSystemConstrain(theSystem,i,value); }}
                           
   return femFullSystemEliminate(theSystem);
}
