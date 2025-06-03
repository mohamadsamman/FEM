clear ; clf
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
//       Icas ::  permet de définir le type de physique. 
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
// =====================================================
function Pk_Bf_at_Xg = Fonctions_de_Base(Nvel, Ngl, Lambda)
    select Nvel 

    case 2 // P1-Lagrange 1D
        // ------------------------------------------------------
        Pk_Bf_at_Xg(1,1:Ngl) =        Lambda(1,1:Ngl)   
        Pk_Bf_at_Xg(2,1:Ngl) =  1.0 - Lambda(1, 1:Ngl); 

    case 3 // P2-Lagrange 1D
        // ------------------------------------------------------

        Pk_Bf_at_Xg(1,1:Ngl) =  Lambda(1,1:Ngl).*( 2.0*Lambda(1,1:Ngl) - 1)   
        Pk_Bf_at_Xg(2,1:Ngl) =  Lambda(2,1:Ngl).*( 2.0*Lambda(2,1:Ngl) - 1)
        Pk_Bf_at_Xg(3,1:Ngl) =  4.0*Lambda(1,1:Ngl).*Lambda(2,1:Ngl)


    case 4 // P3-Lagrange 1D
        // ------------------------------------------------------
        L1 = Lambda(1,1:Ngl)
        L2 = Lambda(2,1:Ngl)
        
        Pk_Bf_at_Xg(1,1:Ngl) =  L1.*( 3.0*L1 - 1).* ( 3.0*L1 - 2)/2  
        Pk_Bf_at_Xg(2,1:Ngl) =  L2.*( 3.0*L2 - 1).* ( 3.0*L2 - 2)/2 
        Pk_Bf_at_Xg(3,1:Ngl) =  9*L1.*L2.*( 3.0*L1 - 1)/2 
        Pk_Bf_at_Xg(4,1:Ngl) =  9*L1.*L2.*( 3.0*L2 - 1)/2 

 // -------------------------------------------------------

    else // P1-Lagrange 1D (par défaut)
        // ------------------------------------------------------
        Pk_Bf_at_Xg(1,1:Ngl) =  Lambda(1,1:Ngl)   
        Pk_Bf_at_Xg(2,1:Ngl) =  1.0 - Lambda(1, 1:Ngl); 

     end

endfunction
// =====================================================
// =====================================================
function Grad_Pk_Bf_at_Xg ...
                  = Grad_de_Fonctions_de_Base(Nvel, Ngl, Lambda, xl1e, xl2e)
    
    Grad_Lambda = MGrad_Lambda_of_P1(xl1e, xl2e) 

    select Nvel 

    case 2 // P1-Lagrange 1D
        // ------------------------------------------------------

        Grad_Pk_Bf_at_Xg(1,1:Ngl) = Grad_Lambda(1)
        Grad_Pk_Bf_at_Xg(2,1:Ngl) = Grad_Lambda(2)


    case 3 // P2-Lagrange 1D
        // ------------------------------------------------------

        Grad_Pk_Bf_at_Xg(1,1:Ngl) =    Grad_Lambda(1)*(4.0*Lambda(1,1:Ngl) - 1)
        Grad_Pk_Bf_at_Xg(2,1:Ngl) =    Grad_Lambda(2)*(4.0*Lambda(2,1:Ngl) - 1) 
        Grad_Pk_Bf_at_Xg(3,1:Ngl) = -4*Grad_Lambda(2)*(2.0*Lambda(2,1:Ngl) - 1)

    case 4 // P3-Lagrange 1D
        // ------------------------------------------------------
        L1 = Lambda(1,1:Ngl)
        L2 = Lambda(2,1:Ngl)
        
        G1 = Grad_Lambda(1)
        G2 = Grad_Lambda(2)

        Grad_Pk_Bf_at_Xg(1,1:Ngl) =    G1.*( 3.0*L1 - 1).*( 3.0*L1 - 2)/2 ...
                                    + L1.*( 3.0*G1    ).*( 3.0*L1 - 2)/2 ...
                                    + L1.*( 3.0*L1 - 1).*( 3.0*G1    )/2;
                                      
        Grad_Pk_Bf_at_Xg(2,1:Ngl) =    G2.*( 3.0*L2 - 1).*( 3.0*L2 - 2)/2 ...
                                    + L2.*( 3.0*G2    ).*( 3.0*L2 - 2)/2 ...
                                    + L2.*( 3.0*L2 - 1).*( 3.0*G2    )/2;
                                      
        Grad_Pk_Bf_at_Xg(3,1:Ngl) =    9*G1.*L2.*( 3.0*L1 - 1)/2 ...
                                    + 9*L1.*G2.*( 3.0*L1 - 1)/2 ...
                                    + 9*L1.*L2.*( 3.0*G1    )/2 

        Grad_Pk_Bf_at_Xg(4,1:Ngl) =    9*G1.*L2.*( 3.0*L2 - 1)/2 ...
                                    + 9*L1.*G2.*( 3.0*L2 - 1)/2 ...
                                    + 9*L1.*L2.*( 3.0*G2    )/2 
// -------------------------------------------------------

    else // P1-Lagrange 1D (par défaut)
        // ------------------------------------------------------
        
        Grad_Pk_Bf_at_Xg(1,1:Ngl) = Grad_Lambda(1)
        Grad_Pk_Bf_at_Xg(2,1:Ngl) = Grad_Lambda(2)

    end

endfunction
// =====================================================
// =====================================================
function [Aex,Bex] = ExactE(Nve, Icas)
    select Nve 
           case 2 // P1-Lagrange 1D
        // ------------------------------------------------------
        // Dans ce contexte, Matrice (Aex) se second membre (Bex) 
        // calculés analytiquement
        // --------------------------------------
        Aex=[ 4/3, -5/6 ;...
        -5/6,  4/3 ]
        Bex=[7/6 ;...
        8/6  ]

    case 3 // P2-Lagrange 1D
        // ------------------------------------------------------
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
        //          
        Aex=[    148. , -13. , -189. ,   54. ;...
                 -13. , 148. ,   54. , -189. ;...
                -189. ,  54. ,  432. , -297.;
                 54.  ,-189. , -297. ,  432.]/40;
        Bex=[17 ; 28 ; 54 ; 81]/120.0;
        // -------------------------------------------------------

    else // P1-Lagrange 1D (par défaut)
        // ------------------------------------------------------

        Aex=[4/3, -5/6 ; -5/6, 4/3]
        Bex=[7/6;8/6]
    end

endfunction

// ==================================================
//--------------------------------------
function out = MyF(icas, x, Mu, Re, Cv, fkl,Ibordl)
   select icas 
    case 1
        out =  x      
    case 2
        out =  x - 1 
    case 10
        a   = - 2*%pi*fkl*cos(2*%pi*fkl*x1g) 
        out =  sin(2*fkl*%pi*x)*( Mu*(2*fkl*%pi)**2 + Re ) 
        if (Ibordl == 2) then
            out =  sin(2*fkl*%pi*x)*( Mu*(2*fkl*%pi)**2 + Re ) + a*Re*(x-3)
        end
        
    else
        out =  sin(2*%pi*x)*( Mu*(2*%pi)**2 + Re )      
    end
endfunction
// =====================================================
// ==================================================
function [Coef_Mu, Coef_Re, Coef_Cv] = Contexte_Physique(Icas)

select Icas
     
case 1  // Réaction-Diffusion
    Coef_Mu = 1.0
    Coef_Re = 1.0
    Coef_Cv = 0.0

case 2 // Diffusion
    Coef_Mu = 1.0
    Coef_Re = 0.0 
    Coef_Cv = 0.0

case 3 // Réaction
    Coef_Mu = 0.0
    Coef_Re = 1.0 
    Coef_Cv = 0.0

else // Par défaut
    Coef_Mu = 1.0
    Coef_Re = 1.0  
    Coef_Cv = 0.0
       
end
endfunction
// ==================================================
// ==================================================
function [Omega,Lambda] = Poids_et_Points_De_Gauss(Ngl)
select  Ngl
case 1 then
    //  donnees pour un point de Gauss
    // -------------------------------------------------
    Omega(1)        = 1.0;
    Lambda(1,1)     = 0.5;
    
case 2 then
    //  donnees pour  deux points de Gauss
    // -------------------------------------------------
    Omega(1:Ngl)    = 0.5;
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

    Lambda(2,1:Ngl) = 1.0 - Lambda(1, 1:Ngl); 

endfunction
// ==================================================
// ==================================================
function out = ExactF(xsl, gsl, Icasl, Ibordl )
    select Icasl 
    case 1
        out = -(xsl - 3).*(xsl.*xsl - gsl )/6
    case 2
        out = -(xsl - 3).*(xsl.*xsl - gsl )/6
    case 3
        out = -(xsl - 3).*(xsl.*xsl - gsl )/6
    case 10 
        if( Ibordl == 1 ) then 
            out = sin(2*fk*%pi*xsl)
        elseif( Ibordl == 2 ) 
            a   = - 2*%pi*fk*cos(2*%pi*fk*x1g)
            out = sin(2*fk*%pi*xsl) + a*(xsl -3)
        end

    end
endfunction

function [CoorSl, NuSl, CoorVl, NuVl] = ...
    Donnees_Interpolation(Nsl, Nel,Nvel, x1gl, x2gl)
    // construction du Maillage géométrique
    // ====================================

    CoorSl       = linspace(x1gl,x2gl,Nsl)

    NuSl(1,1:Nel) = (2:Nsl)'
    NuSl(2,1:Nel) = (1:Nsl-1)'

    // Données de l'interpolation : EF
    // ===============================
    NuVl   = NuSl
    CoorVl = CoorSl

    // ----------------------------------------------------------------------
    //    Données de l'interpolation suivant le type d'élément fini  (Nve)
    // ----------------------------------------------------------------------
    if (Nel == 2) then 
        select Nvel

        case 2
            NuVl   = NuSl
            CoorVl = CoorSl

        case 3
            NuVl   = NuSl
            CoorVl = CoorSl
            is     = Ns

            is  = is + 1
            e   = 1
            l   = 1
            L1  = 1.0/2
            x1l = CoorSl(NuSl(1,e))
            x2l = CoorSl(NuSl(2,e))
            NuVl(2+l,e)  = is
            CoorVl(is)   = L1*x1l + (1-L1)*x2l

            is  = is + 1
            e   = 2
            l   = 1
            L1  = 1.0/2
            x1l = CoorSl(NuSl(1,e))
            x2l = CoorSl(NuSl(2,e))
            NuVl(2+l,e)  = is
            CoorVl(is)   = L1*x1l + (1-L1)*x2l

        else
            /// Cas général pour deux éléments
            // --------------------------------
            NuVl   = NuSl
            CoorVl = CoorSl
            is     = Nsl

            e   = 1
            x1e = CoorSl( NuSl(1,e) ) 
            x2e = CoorSl( NuSl(2,e) ) 
            for l=1:Nvel-2
                is         = is +1
                L1         = 1.0 - l/real(Nvel-1)
                NuVl(l+2,e) = is
                CoorVl(is)  = L1*x1e + (1-L1)*x2e
            end

            e   = 2
            x1e = CoorSl( NuSl(1,e) ) 
            x2e = CoorSl( NuSl(2,e) ) 
            for l=1:Nvel-2
                is         = is +1
                L1         = 1.0 - l/real(Nvel-1)
                NuVl(l+2,e) = is
                CoorVl(is)  = L1*x1e + (1-L1)*x2e
            end
        end

    else
        ///
        //  Cas général 
        // 
        select Nvel

        case 2
            NuVl   = NuSl
            CoorVl = CoorSl
        else
            NuVl   = NuSl
            CoorVl = CoorSl
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            // A compléter
            // ////////////////////
        end
    end 

// ////////////////////////////////////////////////////////////////////////
endfunction
// ///////////////////////////////////////////////////////
//             Fin des fonctions du programme 
// //////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////   
// Définition de l'élément géométrique
// ====================================================
x1g = 2 ; x2g = 3  ; 

// Grad_Lambda = MGrad_Lambda_of_P1(x1e, x2e)

// définition du cas traité
// Par exemple :
//      Icas   = 1   et   Nve = 2    pour "Exercice 2.1.1" EF P1
//      Icas   = 2   et   Nve = 3    pour "Exercice 2.1.2" EF P2
// -------------------------
Icas   = 10   // Physique du problème
fk     = 2
// Icas = 2
//        avec Ibord = 1 : la solution exacte est Uf(x) =(x-3)(x^2 - 4)/6
//        avec Ibord = 2 : la solution exacte est Uf(x) =(x-3)(x^2)/6      
// Icas = 10 
//        avec Ibord = 1 : la solution exacte est Uf(x)= sin(2%pi*fk*x))

Nve    = 3    // type d'élément fini 
//              Nve = 2 ==> pour P1-Lagrange 1D
//              Nve = 3 ==> pour P2-Lagrange 1D
//              Nve = 4 ==> pour P3-Lagrange 1D
//              p refinement quand Nve devient de plus en plus grand 
// Nve = 2,3,4

Ns     = 3   // nombre de points géométriques 
//              h-refinement avec dx de plus en plus petit 
//              ce qui drevient à prendre Ns de plus en plus grand
// Ns =  2, 3, 5, 9, 17 
// ----------------
Ne            = Ns -1           // Nombre d'éléments géométriques
Nv            = Ns + Ne*(Nve-2) // Nombre global de Variables 

// -----------------------------------------------------------

Ibord = 2     // Type de conditions aux limites 
// Ibord = 1 ===>    Dirichlet homogène en x1g et x2g
// Ibord = 2 ===>    Dirichlet homogène en x2g 
//                et Neumann homogène en x1g      
// --------------------------

Ng             = 3 // nombre de points de Gauss pour l'integration numérique.

/////////////////////////////////////////////////////////////////////////////

// -------------------------------------------------------------------
// Définition des points et des poids pour l'intégration Numérique
// -------------------------------------------------------------------
[Omega,Lambda] = Poids_et_Points_De_Gauss(Ng)


// Définition des fonctions de base et gradients aux points de Gauss
// =================================================================
Pk_Bf_at_Xg  = Fonctions_de_Base(Nve, Ng, Lambda)
    

// Définition de la physique du problème 
// ---------------------------------
[Coef_Mu, Coef_Re, Coef_Cv] = Contexte_Physique(Icas)

// ---------------------------------------------------------
// suivant  les paramètres géométriques :Ns, Ne, x1g, x2g
//          et le type d'élément fini : Nve
// construction :
//               des Données géométriques : CoorS, NuS
//               des Données d'Interpolation 
// 
// ----------------------------------------------------------------------
[CoorS, NuS, CoorV, NuV] = Donnees_Interpolation(Ns, Ne, Nve, x1g, x2g)


// ////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////
// Calcul global de la Matrice (Mat_A) et du second membre (Vec_B)
// //////////////////////////////////////////////////////////////////////////-----------------------
// Initialisations à zéro
Mat_A(1:Nv,1:Nv)  = 0
Vec_B(1:Nv)       = 0
//////////////////////////////////////////////////////////////////////
//                   Assemblage par éléménts 
// -------------------------------------------------------------------
// Boucle sur les éléments 
//---------------------------------
for e = 1:Ne 
    // données géométriques locales à l'élément.
    // -----------------------------------------
    x1e = CoorS( NuS(1,e) ) 
    x2e = CoorS( NuS(2,e) ) 
    // Jacobien de la transformation entre (x1e, x2e) et [0,1]
    JacE         = abs(x2e-x1e);
        
    Grad_Pk_Bf_at_Xg = Grad_de_Fonctions_de_Base(Nve, Ng, Lambda, x1e, x2e)

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
            is        = NuV(l,e)
            Vl        =      Pk_Bf_at_Xg(l,g)
            GVl       = Grad_Pk_Bf_at_Xg(l,g)

            fg        = MyF(Icas, xg, Coef_Mu, Coef_Re, Coef_Cv,fk,Ibord)

            bl        = JacE*fg*Vl 

            Vec_B(is) =  Vec_B(is)  + Omega(g)*bl 


            // seconde boucle (inbriquée) sur les ddl locaux
            // =============================================
            // Assemblage local de la matrice 
            for m = 1:Nve
                js       = NuV(m,e)
                Vm       =      Pk_Bf_at_Xg(m,g)
                GVm      = Grad_Pk_Bf_at_Xg(m,g)

                Alm       = JacE * ( Coef_Mu*GVl*GVm + Coef_Re*Vl*Vm )

                Mat_A(is,js) =  Mat_A(is,js) + Omega(g)*Alm 
                
            end // Fin de la boucle sur m
            
        end // Fin de la boucle sur l
        
    end // Fin de la boucle sur les points de Gauss
    
end // Fin de la boucle sur les éléments
// --------------------------------------
// ////////////////////////////////////////////////////////////////////////

// Matrice et second membre calculés analytiquement
// -------------------------------------------------
[Aex,Bex] = ExactE(Nve, Icas)

// Affichages 
 disp("Mat_A", Mat_A)
 disp("Vec_B", Vec_B)
//disp("norm(Aex - Mat_Ae)", norm(Aex - Mat_A))
//disp("norm(Bex - Vec_Be)", norm(Bex - Vec_B))
// ////////////////////////////////////////////////////////////////////////
//                      Conditions aux limites 
// ////////////////////////////////////////////////////////////////////////
// -----------------------
//     Ibord = 1 ===>    Dirichlet homogène en x1e et x2e
//     Ibord = 2 ===>    Dirichlet homogène en x2e 
//                    et Neumann homogène en x1e
// -----------------------------------------------------------
//////////
select Ibord
    
case 1 // Conditions aux limites de Dirichlet homogène 
    is = NuS(2,1)
    Mat_A(is,:) = 0 ; Mat_A(is,is) = 1 ; Vec_B(is) = 0
    
    is = NuS(1,Ne)
    Mat_A(is,:) = 0 ; Mat_A(is,is) = 1 ; Vec_B(is) = 0  
    
    ExactU=[0;0;3/16]
    gm    = 4

case 2 // Conditions de Dirichlet (l=2) et Neumann (l=1) homogène 
    is = NuS(1,Ne)
    Mat_A(is,:) = 0 ; Mat_A(is,is) = 1 ; Vec_B(is) = 0  
    
    ExactU=[2/3;0;25/48]
    gm = 0
 
case 3 // autres conditions aux limites à définir
    
else  // Conditions aux limites par défaut
 Mat_A(1,:) = 0 ; Mat_A(1,1) = 1 ; Vec_B(1) = 0
 Mat_A(2,:) = 0 ; Mat_A(2,2) = 1 ; Vec_B(2) = 0 
      
end
// //////////////////////////////////////////////////////////////////////
//         Solution (Uh) du système algébrique (linéaire) 
// //////////////////////////////////////////////////////////////////////
//
Uh = Mat_A\Vec_B

// //////////////////////////////////////////////////////////////////////
//         Calcul des erreurs entre la fonction exacte (Uf)
//         et la fontion interpolée approchée  (Uhf) 
// //////////////////////////////////////////////////////////////////////
//
rf  = 1
ErrL1(rf) = 0;   ErrL2(rf)   = 0; ErrLinf(rf) = 0.0 ;

for e = 1:Ne 
    x1e = CoorS( NuS(1,e) ) 
    x2e = CoorS( NuS(2,e) ) 
    // Jacobien de la transformation entre (x1e, x2e) et [0,1]
    JacE         = abs(x2e-x1e);
    // Boucle sur les points de Gauss
    // ----------------------------------------------------------
    for g = 1:Ng
        // Fonctions de base Pk
        // *****************************************
        // Phi  = Phi_Pk(Lambda_g(1,gp), Lambda_g(2,gp), EF_Pk)
        wg   = Omega(g);

        // Point physique associé au point de Gauss
        // ------------------------------------------        
        xg  = Lambda(1,g)*x1e + Lambda(2,g)*x2e ;

        // Solution approchée au point de Gauss
        // ------------------------------------------        
        Uhg = 0
        for l=1:Nve 
            is    = NuV(l,e)  
            Vl    =   Pk_Bf_at_Xg(l,g)           
            Uhg   = Uhg +   Uh(is)*Vl;
        end 

        // Solution exacte au point de Gauss
        // -----------------------------------------
        Uex  = ExactF(xg,gm,Icas,Ibord)    ;

        // Les Erreurs L1, L2 et Linf
        // ----------------------------------------------
        ErrL1(rf)   = ErrL1(rf) + wg*JacE*abs(Uhg -Uex)  ; 
        ErrL2(rf)   = ErrL2(rf) + wg*JacE*(Uhg -Uex)^2  ;
        ErrLinf(rf) = max( ErrLinf(rf), abs(Uhg -Uex) ) ; 
    end 
end 
  
ErrL2(rf) = sqrt(ErrL2(rf))

// Affichage des erreurs dans la fenetre console
// -----------------------------------------------
disp("Erreur en")
disp( ["Norme L1", "Norme L2", "Norme Linf"] )
disp([ ErrL1; ErrL2; ErrLinf' ]')


// -------------------------------------------------------------
// Pour tracer la solution approchée et la solution exacte
// -------------------------------------------------------------
Nplot      = 100
Xplot(1,:) = linspace(0,1,Nplot)
Xplot(2,:) = 1 - Xplot(1,:)
//
// fonctions de base et gradients aux points Xplot
// -------------------------------------------------------------------------
Plot_Bf = Fonctions_de_Base(Nve, Nplot, Xplot )
// -------------------------------------------------------------------------

// Fonctions approchée (Uhf)  et exact (Uf) aux points Xplot
// ---------------------------------------------------------
Uhf = zeros(1, Ne*Nplot)
for e=1:Ne
    x1e = CoorS( NuS(1,e) ) 
    x2e = CoorS( NuS(2,e) )

    // Fonction exacte (Uf) aux points Xplot
// ------------------------------------------
xs  = x1e*Xplot(1,:) + x2e*Xplot(2,:)

is1 = 1   + (e-1)*Nplot
is2 = is1 - 1 + Nplot

Uf(is1:is2)  = ExactF(xs,gm,Icas,Ibord)
Xs(is1:is2)  = xs

for l=1:Nve 
    is = NuV(l,e)          
    Uhf(is1:is2)  = Uhf(is1:is2) +  Uh(is)*Plot_Bf(l,:)
end 
end
// ------------------------------------------

/////////////////////////////////////////////////////////////////////////////
// Ici on trace la solution exacte (Uf) et la solution approchée (Uhf)
// -------------------------------------------------------------------
plot(Xs   , Uf      ,"-r" ,'LineWidth',3 )
plot(Xs   , Uhf     ,"-b",'LineWidth',2 )

// Solution approchée aux points d'interpolation (en magenta)
// ----------------------------------------------------------
plot(CoorV, Uh'      , "om",'LineWidth',2 )

// Solution approchée aux points géométriques (en bleu)
// -----------------------------------------------------
plot(CoorS, Uh(1:Ns)', "*b",'LineWidth',3 )
//------------------------------------------

//////////////////////////////////////////////////////////////////////////////
