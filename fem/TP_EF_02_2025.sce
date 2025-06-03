clear 
// ///////////////////////////////////////////////////////////////////////////
//  Programme de calcul numérique sur un élément défini par x1, x2
//  De la matrice et du second membre locaux à l'élément. 
//   Ng,  Lambda(1:Ng, 1) == Lambda_1  et Omega(1:Ng) 
//   on en déduit 
//    Lambda(1:Ng, 2)  = 1 -  Lambda(1:Ng, 1) == Lambda_2
// ----------------------------------------------------------------
//  les données du programe sont:
//       x1, x2
//       Ng, omega(1:Ng), Lambda1(1:Ng) et on en déduit  Lambda2(1:Ng)
//       1) La fonction évaluant les fonctions de base P1 
//       2) La fonction évaluant les gradients des fonctions de base P1
//       3) La fonction definissant le second membre du problème 
//       Nve  ::  permet de choisir le type d'élément fini. 
//       Icas ::  permet de définir le type de second membre. 
//  --------------------------------------------------------------------
// ///////////////////////////////////////////////////////////////////////////

// ///////////////////////////////////////////////////////
//             Les fonctions du programme 
// ///////////////////////////////////////////////////////

// --------------------------------------------------------------
// fonctions de base P1 
//-------------------------------------------------------------
function out = Lambda_of_P1(x, x1, x2)
   out(1) = (x - x2)/(x1-x2)
   out(2) = (x - x1)/(x2-x1)
endfunction

//--------------------------------
// Gradients des fonctions de Base P1
//--------------------------------
function G = MGrad_Lambda_of_P1(x1, x2)  
    G(1) = 1.0/(x1-x2)
    G(2) = 1.0/(x2-x1)
endfunction
//---------------------------------------------------------------------
//   fonction utilisée pour la validation de l'intégration Numérique
// -------------------------------------------------------------------- 
//--------------------------------------
// fonction second membre de problème.
//--------------------------------------
function out = MyF(icas, x)
   select icas 
    case 1
        out =  x      
    case 2
        out =  x - 1  
    else
        out =  x      
    end
endfunction

// ///////////////////////////////////////////////////////
//             Fin des fonctions du programme 
// ///////////////////////////////////////////////////////

// -------------------------------------------------------------------
// Définition des points et des poids pour l'intégration Numérique
// -------------------------------------------------------------------
Ng    = 3 // nombre de points de Gauss pour l'integration numérique.
select  Ng
case 1 then
    //  donnees pour un point de Gauss
    // -------------------------------------------------
    Omega(1)        = 1.0;
    Lambda(1,1)     = 0.5;
    
case 2 then
    //  donnees pour  deux points de Gauss
    // -------------------------------------------------
    Omega(1:Ng)    = 0.5;
    Lambda(1,1)    = 0.5*( 1 - 1.0/sqrt(3)) ;
    Lambda(1,2)    = 0.5*( 1 + 1.0/sqrt(3)) ;
    
case 3 then
    //  donnees pour trois points de Gauss
    // ------------------------------------------------- 

    Omega(1)  = 5.0/18.0 ;
    Omega(2)  = 8.0/18.0 ;
    Omega(3)  = 5.0/18.0 ;

    Lambda(1, 1) = 0.5 * ( 1.0 - sqrt(3.0/5.0) ) ;
    Lambda(1, 2) = 0.5   ;
    Lambda(1, 3) = 0.5 * ( 1.0 + sqrt(3.0/5.0) )   ;

end 

    Lambda(2,1:Ng) = 1.0 - Lambda(1, 1:Ng); 
// ////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////   
// Définition de l'élément géométrique
// ====================================================
x1e = 2 ; x2e = 3  ; 

// définition du cas traité
// Par exemple :
//      Icas   = 1   et   Nve = 2    pour "Exercice 2.1.1" EF P1
//      Icas   = 2   et   Nve = 3    pour "Exercice 2.1.2" EF P2
// -------------------------
Icas   = 2   // Physique du problème
Nve    = 3   // type d'élément fini  
// --------------------------

Grad_Lambda = MGrad_Lambda_of_P1(x1e, x2e)

// fonctions de base et gradients aux points de Gauss
// Nve = 2 ==> pour P1-Lagrange 1D
// Nve = 3 ==> pour P2-Lagrange 1D
// Nve = 4 ==> pour P3-Lagrange 1D
// ====================================================

select Nve 
    
case 2 // P1-Lagrange 1D
    // ------------------------------------------------------
    Pk_Bf_at_Xg(1,1:Ng) =        Lambda(1,1:Ng)   
    Pk_Bf_at_Xg(2,1:Ng) =  1.0 - Lambda(1, 1:Ng); 

    Grad_Lambda = MGrad_Lambda_of_P1(x1e, x2e) 

    Grad_Pk_Bf_at_Xg(1,1:Ng) = Grad_Lambda(1)
    Grad_Pk_Bf_at_Xg(2,1:Ng) = Grad_Lambda(2)

    // Dans ce contexte, Matrice (Aex) se second membre (Bex) 
    // calculés analytiquement
    // --------------------------------------
    Aex=[ 4/3, -5/6 ;...
         -5/6,  4/3 ]
    Bex=[7/6 ;...
         8/6  ]

case 3 // P2-Lagrange 1D
    // ------------------------------------------------------

    Pk_Bf_at_Xg(1,1:Ng) =  Lambda(1,1:Ng).*( 2.0*Lambda(1,1:Ng) - 1)   
    Pk_Bf_at_Xg(2,1:Ng) =  Lambda(2,1:Ng).*( 2.0*Lambda(2,1:Ng) - 1)
    Pk_Bf_at_Xg(3,1:Ng) =  4.0*Lambda(1,1:Ng).*Lambda(2,1:Ng)

    Grad_Lambda = MGrad_Lambda_of_P1(x1e, x2e) 

    Grad_Pk_Bf_at_Xg(1,1:Ng) =    Grad_Lambda(1)*(4.0*Lambda(1,1:Ng) - 1)
    Grad_Pk_Bf_at_Xg(2,1:Ng) =    Grad_Lambda(2)*(4.0*Lambda(2,1:Ng) - 1) 
    Grad_Pk_Bf_at_Xg(3,1:Ng) = -4*Grad_Lambda(2)*(2.0*Lambda(2,1:Ng) - 1)
           
    // Dans ce contexte, Matrice (Aex) se second membre (Bex) 
    // calculés analytiquement : pour 
    // --------------------------------------
    Aex=[ 7/3,  1/3, -8/3 ; ...
          1/3,  7/3, -8/3 ; ...
         -8/3, -8/3, 16/3  ];
    Bex=[1/6 ; ...
         2/6 ; ...
         1];

case 4 // P3-Lagrange 1D
    // ------------------------------------------------------
    //          À faire
    // -------------------------------------------------------

else // P1-Lagrange 1D (par défaut)
    // ------------------------------------------------------
    Pk_Bf_at_Xg(1,1:Ng) =  Lambda(1,1:Ng)   
    Pk_Bf_at_Xg(2,1:Ng) =  1.0 - Lambda(1, 1:Ng); 

    Grad_Lambda = MGrad_Lambda_of_P1(x1e, x2e) 

    Grad_Pk_Bf_at_Xg(1,1:Ng) = Grad_Lambda(1)
    Grad_Pk_Bf_at_Xg(2,1:Ng) = Grad_Lambda(2)

    Aex=[4/3, -5/6 ; -5/6, 4/3]
    Bex=[7/6;8/6]
end

// choix de la physique du problème 
// ---------------------------------
select Icas
     
case 1  // Réaction-Diffusion
    Coef_Mu = 1.0
    Coef_Re = 1.0

case 2 // Diffusion
    Coef_Mu = 1.0
    Coef_Re = 0.0 
    
case 3 // Réaction
    Coef_Mu = 0.0
    Coef_Re = 1.0 
    
else // Par défaut
    Coef_Mu = 1.0
    Coef_Re = 1.0  
           
end
// Calcul local de la Matrice (Mat_Ae) et du second membre (Vec_Be)
// ------------------------------------------------------------------

// Initialisations à zéro
//-----------------------
Mat_Ae(1:Nve,1:Nve) = 0
Vec_Be(1:Nve)       = 0

// Jacobien de la transformation entre (x1e, x2e) et [0,1]
JacE         = abs(x2e-x1e);

// boucles sur les points de Gauss
// ===============================
for g = 1:Ng

    // Image du point de Gauss dans xsi_g dans (x1e, x2e)
    // --------------------------------------------------
    xg = Lambda(1,g)*x1e + Lambda(2,g)*x2e ;
    
    // -----------------------------------------------------
    // Assemblage local du second membre et de la matrice 
    // -------------------------------------------------------
    
    // première boucle sur le ddl locaux
    // Assemblage local du second membre
    // ==================================
    for l = 1:Nve

        Vl        =      Pk_Bf_at_Xg(l,g)
        GVl       = Grad_Pk_Bf_at_Xg(l,g)

        fg        = MyF(Icas, xg)

        bl        = JacE*fg*Vl 

        Vec_Be(l) =  Vec_Be(l)  + Omega(g)*bl 
        

        // seconde boucle (inbriquée) sur les ddl locaux
        // =============================================
        // Assemblage local de la matrice 
        for m = 1:Nve
            Vm       =      Pk_Bf_at_Xg(m,g)
            GVm      = Grad_Pk_Bf_at_Xg(m,g)
            
            Alm       = JacE * ( Coef_Mu*GVl*GVm + Coef_Re*Vl*Vm )

            Mat_Ae(l,m) =  Mat_Ae(l,m) + Omega(g)*Alm 
        end
    end
end
disp(norm(Aex - Mat_Ae))
disp(norm(Bex - Vec_Be))
// ------------------------------------------------------------
// Conditions aux limites 
// -----------------------
//     Ibord = 1 ===>    Dirichlet homogène en x1e et x2e
//     Ibord = 2 ===>    Dirichlet homogène en x2e 
//                    et Neumann homogène en x1e
// -----------------------------------------------------------

/*
Ibord = 2
select Ibord
    
case 1 // Conditions aux limites de Dirichlet homogène 
 Mat_Ae(1,:) = 0 ; Mat_Ae(1,1) = 1 ; Vec_Be(1) = 0
 Mat_Ae(2,:) = 0 ; Mat_Ae(2,2) = 1 ; Vec_Be(2) = 0  
 Exact=[0;0;3/16]
 
case 2 // Conditions de Dirichlet (l=2) et Neumann (l=1) homogène 
 Mat_Ae(2,:) = 0 ; Mat_Ae(2,2) = 1 ; Vec_Be(2) = 0
 Exact=[2/3;0;25/48]

case 3 // autres conditions aux limites à définir
    
else  // Conditions aux limites par défaut
 Mat_Ae(1,:) = 0 ; Mat_Ae(1,1) = 1 ; Vec_Be(1) = 0
 Mat_Ae(2,:) = 0 ; Mat_Ae(2,2) = 1 ; Vec_Be(2) = 0 
      
end

Resu = Mat_Ae\Vec_Be
// disp(Resu)
disp(norm(Resu-Exact))
*/

