clear
clf

function [out] = Exact(x, y, gamx, gamy)
    out = sin(gamx*x).*cos(gamy*y)
endfunction

function [out] = f(Xg, gamx, gamy, betx, bety, cx, cy, alpha)
    x = Xg(1)
    y = Xg(2)

    term1 = (alpha + betx * gamx^2 + bety * gamy^2) * sin(gamx * x) * cos(gamy * y);
    term2 = cx * gamx * cos(gamx * x) * cos(gamy * y);
    term3 = -cy * gamy * sin(gamx * x) * sin(gamy * y);

    out = term1 + term2 + term3; 
 endfunction
//--------------------------------
function [betx, bety , cx, cy, alpha] = Composite_Mat(Xg,choix)
    // Sélection des paramètres selon le choix
    select choix
    case 1 then
        // *** Problème de réaction ***
        // betx = 0, bety = 0, cx = 0, cy = 0, alpha ≠ 0
        betx = 0.0
        bety = 0.0
        cx    = 0.0
        cy    = 0.0
        alpha = 1.1057110571

    case 2 then
        // *** Problème de diffusion isotrope ***
        // betx = bety > 0, cx = 0, cy = 0, alpha = 0
        betx = 1.1057110571  // Exemple de valeur
        bety = 1.1057110571  // betx = bety
        cx    = 0.0
        cy    = 0.0
        alpha = 0.0

    case 3 then
        // *** Problème de diffusion anisotrope ***
        
        // betx > 0, betx ≠ bety > 0, cx = 0, cy = 0, alpha = 0
        betx = 2.1057110571  // Exemple de valeur
        bety = 1.1057110571  // betx ≠ bety
        
        // betx > 0,  bety = 0, cx = 0, cy = 0, alpha = 0
        //betx = 2.1057110571
        //bety = 0.0
        
        // betx = 10e-8, betx = 1, cx = 0, cy = 0, alpha = 0
        //betx = 0.00000008  
        //bety = 1.0

        cx    = 0.0
        cy    = 0.0
        alpha = 0.0

    case 4 then
        // *** Problème de convection ***
        // betx = 0, bety = 0, cx ≠ 0, cy ≠ 0, alpha = 0
        betx = 0.0
        bety = 0.0
        cx    = 1.1057110571  // Exemple de valeur
        cy    = 1.1057110571  // Exemple de valeur
        alpha = 0.0

    case 5 then
        // *** Problème de convection-réaction-diffusion anisotrope 2D ***
        // betx = 1, bety = 2, cx = 1, cy = 0.5, alpha = -5
        betx = 1.0
        bety = 2.0
        cx    = 1.0
        cy    = 0.5
        alpha = -5.0

    else
        // Choix invalide : valeurs par défaut
        disp("Choix non reconnu. Utilisation des paramètres par défaut.")
        betx = 0.0
        bety = 0.0
        cx    = 0.0
        cy    = 0.0
        alpha = 1.0
    end
endfunction
//--------------------------------
// Construction d'un maillage Pk de triangles. 
//--------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------
function [N_vertives, N_elements, Coor, Nu, LogP, LogE, NuV] ...
                          = Struct_Pk_Mesh(Nx, Ny, Lx, Ly,pk)
    
    MNx        = Nx + (pk-1)*(Nx - 1) 
    MNy        = Ny + (pk-1)*(Ny - 1) 
    
    N_vertives = MNx*MNy
    N_elements = 2*(Nx-1)*(Ny-1)
    
    dx = Lx/(MNx -1)
    dy = Ly/(MNy -1)
    is = 0
    
    for i=1:MNx
        xi = (i-1)*dx
        for j=1:MNy
            yj = (j-1)*dy
            is = is +1
            Coor(1:2, is)=[xi,yj]'
            LogP(is) = 0
            
            if i==1 then
                LogP(is) = -1
            elseif i==MNx then
                LogP(is) = -2
            elseif j==1 then
                LogP(is) = -10
            elseif j==MNy then
                LogP(is) = -20
            end
        end
    end
    
   ie = 0
    for i=1:Nx-1
        for j=1:Ny-1
            
            is = 1 + pk*(j-1) + (i-1)*pk*MNy ; 
            js = is + pk*MNy ;
            ks = js + pk ;
            ps = is + pk ;
            
            // if( 0 == 0 ) then
                // element du bas 
                // is ---> js  ---> ps ---> is
                // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ie = ie +1
                
                Nu(1:3, ie) = [is  ; js ;  ps] ;
                
                il  = 3 ; 
                iie = 2
                // horizontale  (is --> js)
                // suivant le vecteur (1,0)
                // -----------------------
                for lx=1:pk-1
                    il         = il +1
                    Nu(il, ie) = is + lx*MNy ; 
                end 
                iie = il - pk + 2
                ije = il
                // oblique   (js --> ps ) 
                // suivant le vecteur (-1,1)
                // ----------------------------------
                for lo=1:pk-1
                    il         = il +1
                    Nu(il, ie) = js  + lo - lo*MNy ; 
                end 
                jje= il - pk + 2
                jke=il
                // verticale   (ps --> is)
                // suivant le vecteur (0,-1)
                // -----------------------------
                for ly=1:pk-1
                    il         = il +1
                    Nu(il, ie) = ps - ly ; 
                end 
                kke=il - pk +2
                kie=il
                // -----------------------------
                // ddl dans l'élément du bas
                // ----------------------------
                // -----------------------------
                for lx =1:pk-2
                    for ly = 1:lx
                        il = il +1
                        Nu(il, ie) = is + lx*MNy + ly
                    end
                end
                
               
            
                LogE(ie)    = 0
                            
                // element du haut 
                // js ---> ks ---> ps ---> js
                // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ie = ie +1 ;
                
                Nu(1:3, ie) = [js ;  ks ; ps] ;
                
                il = 3 ; 
                // verticale  (js --> ks)
                // suivant le vecteur (0,+1)
                // -----------------------------
                for ly=1:pk-1
                    il         = il +1
                    Nu(il, ie) = js + ly ; 
                end 
                // horizontale  (ks --> ps)
                // suivant le vecteur (0,-1)
                // -----------------------------
                for lx=1:pk-1
                    il         = il +1
                    Nu(il, ie) = ks  - lx*MNy ; 
                end 
                // oblique   (ps --> js ) 
                // suivant le vecteur (-1,1)
                // ----------------------------------
                for lo=1:pk-1
                    il         = il +1
                    Nu(il, ie) = ps - lo + lo*MNy; 
                end 
                 
                // -----------------------------
                // ddl dans l'élément du haut
                // ----------------------------
                // -----------------------------
                for ly =1:pk-2
                    for lx = 1:ly
                        il = il +1
                        Nu(il, ie) = js - lx*MNy + ly +1
                    end
                end
                
                LogE(ie)    = 0
        
        end
    end
    
     
                    
   
    if (pk==1) then
        NuV(1:3,1) = [1;2;3]
        
    elseif (pk==2) then
        NuV(1:3,1) = [1;4;6]
        NuV(1:3,2) = [4;5;6]
        NuV(1:3,3) = [4;2;5]
        NuV(1:3,4) = [6;5;3]
        
    elseif (pk==3) then
        NuV(1:3,1) = [1 ; 4; 9]
        NuV(1:3,2) = [4 ;10; 9]
        NuV(1:3,3) = [4 ; 5;10]
        NuV(1:3,4) = [5 ; 6;10]
        NuV(1:3,5) = [5 ; 2; 6]
        NuV(1:3,6) = [9 ;10; 8]
        NuV(1:3,7) = [10; 7; 8]
        NuV(1:3,8) = [10; 6; 7]
        NuV(1:3,9) = [8 ; 7; 3]
     end
    

endfunction
//---------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------
//
// Fonctions de base Pk-Lagrange pour k= 1, 2 et 3
// --------------------------------------------
function Phi = Phi_Pk( L1,   L2, pk)
    L3 = 1.0 -   L1  -   L2
    select pk

    case 1 then
        Phi(1) =   L1 ; 
        Phi(2) =   L2 ;
        Phi(3) =   L3 ;

    case 2 then
        Phi(1) =   L1.*( 2.0*  L1 - 1.0) ; 
        Phi(2) =   L2.*( 2.0*  L2 - 1.0) ;
        Phi(3) =   L3.*( 2.0*  L3 - 1.0) ;    
        Phi(4) = 4.0*  L1.*  L2
        Phi(5) = 4.0*  L2.*  L3
        Phi(6) = 4.0*  L3.*  L1

    case 3 then
        Phi(1) = 0.5 * L1 .* (3*L1 - 1) .* (3*L1 - 2); 
        Phi(2) = 0.5 * L2 .* (3*L2 - 1) .* (3*L2 - 2);
        Phi(3) = 0.5 * L3 .* (3*L3 - 1) .* (3*L3 - 2);
        Phi(4) = (9/2) * L1 .* L2 .* (3*L1 - 1);
        Phi(5) = (9/2) * L1 .* L2 .* (3*L2 - 1);
        Phi(6) = (9/2) * L3 .* L2 .* (3*L2 - 1);
        Phi(7) = (9/2) * L2 .* L3 .* (3*L3 - 1);
        Phi(8) = (9/2) * L3 .* L1 .* (3*L3 - 1);
        Phi(9) = (9/2) * L3 .* L1 .* (3*L1 - 1);
        Phi(10) = 27 * L1 .* L2 .* L3;
    else
        disp("Test functions not yet available for EF_Pk = ", pk);
        abort

    end
endfunction
//
// Gradients des fonctions de base Pk-Lagrange
// pour k= 1, 2 et 3
// --------------------------------------------
function GradPhi = GradPhi_Pk(  L1,   L2, G1, G2, pk)
    L3 = 1.0 -   L1  -   L2;
    G3   = - G1 - G2;
    select pk
    case 1 then
        GradPhi(1,:) = G1 ; 
        GradPhi(2,:) = G2 ;
        GradPhi(3,:) = G3 ;
    case 2 then
        GradPhi(1,:) =     G1.*( 2.0*  L1 - 1.0) ...
                       +   L1.*( 2.0*G1        ) ; 
        GradPhi(2,:) =     G2.*( 2.0*  L2 - 1.0) ...
                       +   L2.*( 2.0*G2        ) ;
        GradPhi(3,:) =     G3.*( 2.0*  L3 - 1.0) ...
                       +   L3.*( 2.0*G3        ) ;
        GradPhi(4,:) =  4.0*( G1.*  L2  +   L1.*G2 )
        GradPhi(5,:) =  4.0*( G2.*  L3  +   L2.*G3 )
        GradPhi(6,:) =  4.0*( G3.*  L1  +   L3.*G1 )

    case 3 then
        // φ₁ = 0.5 * L1 * (3*L1 - 1) * (3*L1 - 2)
        GradPhi(1,:) = 0.5 * G1 .* (3*L1 - 1) .* (3*L1 - 2) + 0.5 * L1 .* 3*G1 .* (3*L1 - 2) + 0.5 * L1 .* (3*L1 - 1) .* 3*G1;

        // φ₂ = 0.5 * L2 * (3*L2 - 1) * (3*L2 - 2)
        GradPhi(2,:) = 0.5 * G2 .* (3*L2 - 1) .* (3*L2 - 2) + 0.5 * L2 .* 3*G2 .* (3*L2 - 2) + 0.5 * L2 .* (3*L2 - 1) .* 3*G2;

        // φ₃ = 0.5 * L3 * (3*L3 - 1) * (3*L3 - 2)
        GradPhi(3,:) = 0.5 * G3 .* (3*L3 - 1) .* (3*L3 - 2) + 0.5 * L3 .* 3*G3 .* (3*L3 - 2) + 0.5 * L3 .* (3*L3 - 1) .* 3*G3;

        // φ₄ = (9/2) * L1 * L2 * (3*L1 - 1)
        GradPhi(4,:) = (9/2) * (G1 .* L2 .* (3*L1 - 1) + L1 .* G2 .* (3*L1 - 1) + L1 .* L2 .* 3*G1);

        // φ₅ = (9/2) * L1 * L2 * (3*L2 - 1)
        GradPhi(5,:) = (9/2) * (G1 .* L2 .* (3*L2 - 1) + L1 .* G2 .* (3*L2 - 1) + L1 .* L2 .* 3*G2);

        // φ₆ = (9/2) * L2 * L3 * (3*L2 - 1)
        GradPhi(6,:) = (9/2) * (G2 .* L3 .* (3*L2 - 1) + L2 .* G3 .* (3*L2 - 1) + L2 .* L3 .* 3*G2);

        // φ₇ = (9/2) * L2 * L3 * (3*L3 - 1)
        GradPhi(7,:) = (9/2) * (G2 .* L3 .* (3*L3 - 1) + L2 .* G3 .* (3*L3 - 1) + L2 .* L3 .* 3*G3);

        // φ₈ = (9/2) * L3 * L1 * (3*L3 - 1)
        GradPhi(8,:) = (9/2) * (G3 .* L1 .* (3*L3 - 1) + L3 .* G1 .* (3*L3 - 1) + L3 .* L1 .* 3*G3);

        // φ₉ = (9/2) * L3 * L1 * (3*L1 - 1)
        GradPhi(9,:) = (9/2) * (G3 .* L1 .* (3*L1 - 1) + L3 .* G1 .* (3*L1 - 1) + L3 .* L1 .* 3*G1);

        // φ₁₀ = 27 * L1 * L2 * L3
        GradPhi(10,:) = 27 * (G1 .* L2 .* L3 + L1 .* G2 .* L3 + L1 .* L2 .* G3);

    else

        disp("not yet available for EF_Pk = ", pk)
    end
endfunction
//-----------------------------------------
// déterminant de deux vecteurs 2D
//-----------------------------------------
function out = Determinant(p,q)
    out = p(1)*q(2) - p(2)*q(1)
endfunction

//----------------------------------------- 
// Coordonnées barycentriques 
// au point Xc  dans le triangle ( Xi, Xj, Xk )
// -----------------------------------------------------------------
function Lambda=Lambda_of_P1(Xc, Xi, Xj, Xk)
    Lambda(1) = Determinant(Xc-Xj,Xk-Xj)/Determinant(Xi-Xj,Xk-Xj);
    Lambda(2) = Determinant(Xc-Xk,Xi-Xk)/Determinant(Xj-Xk,Xi-Xk);
    Lambda(3) = Determinant(Xc-Xi,Xj-Xi)/Determinant(Xk-Xi,Xj-Xi);
endfunction
//------------------------------------------------------------------
// Gradient des Coordonnées barycentriques
//--------------------------------------------------------------
function [Grad]=MGrad_Lambda_of_P1(Xi, Xj, Xk)
    
     e1=[1,0] ; e2 =[0,1]
    
    Grad(1,1) = Determinant(e1,Xk-Xj)/Determinant(Xi-Xj,Xk-Xj);
    Grad(1,2) = Determinant(e2,Xk-Xj)/Determinant(Xi-Xj,Xk-Xj);  

    Grad(2,1) = Determinant(e1,Xi-Xk)/Determinant(Xj-Xk,Xi-Xk);
    Grad(2,2) = Determinant(e2,Xi-Xk)/Determinant(Xj-Xk,Xi-Xk);

    Grad(3,1) = Determinant(e1,Xj-Xi)/Determinant(Xk-Xi,Xj-Xi);
    Grad(3,2) = Determinant(e2,Xj-Xi)/Determinant(Xk-Xi,Xj-Xi);
  
endfunction
//---------------------------------------------------------------
//--------------------------------
//--------------------------------

//---------------------------------------------------------------
//---------------------------------------------------------------
// points de Gauss et poids de Gauss 
//   dans l'élément de référence
//---------------------------------------------------------------
//---------------------------------------------------------------
//------------------------------------------------------------
// six points de Gauss
// -------------------------------------------------------
// Nombre de points de Gauss
// --------------------------
Ng                 = 6;
 
 s1 = 0.11169079483905;
 s2 = 0.0549758718227661;
 aa = 0.445948490915965;
 bb = 0.091576213509771;
// poids de Gauss 
// -------------------------------------
 W(1:3)             = s2; 
 W(4:6)             = s1;

// Points de Gauss xi_1 == Lambda_1
// -------------------------------------
 Lambda(1,1:Ng)         = [bb, 1-2*bb, bb, aa, aa, 1-2*aa];
// Points de Gauss xi_1 == Lambda_2
// -------------------------------------
 Lambda(2,1:Ng)         = [bb, bb, 1-2*bb, 1-2*aa, aa, aa];

// Lambda_3 = 1 - Lambda_1 -Lambda_2
// ---------------------------------------
Lambda(3,1:Ng)= 1.0 - Lambda(1,1:Ng) - Lambda(2,1:Ng)


//
//=========================================================
//   Choix de la méthode Eléments finis : 2D, Pk-Lagrange
//--------------------------------------------------------
//   Pk  avec k = EF_Pk 
// ------------------------
//   EF_Pk = 1 <===> P1-Lagrange
//   EF_Pk = 2 <===> P2-Lagrange
//   EF_Pk = 3 <===> P3-Lagrange
//=========================================================
EF_Pk     = 3;
//=========================================================
//   Choix du cas 
//--------------------------------------------------------
//   choix = 1 <===> Problème de réaction
//   choix = 2 <===> Problème de diffusion isotrope
//   choix = 3 <===> Problème de diffusion anisotrope
//   choix = 4 <===> Problème de convection
//   choix = 5 <===> Problème de convection-réaction-diffusion anisotrope 2D
//=========================================================
choix = 1;


// définition des paramètres du maillage.
// ---------------------------------------


Lx  =1   ; Ly=1 

gamx = 1.0*%pi/Lx
gamy = 2.0*%pi/Ly



// Analyse de convergence
// -----------------------------
Nconv    = 4
Vconv    = [13,25, 31, 61, 73, 91; ...
             7,13, 16, 31, 37, 46; ...
             5, 9, 11, 21, 25, 31  ]
//[5,11,21,31,9,21,31,61,17,31,61,121]
for rf = 1:Nconv

    sf =rf
    Nx = Vconv(EF_Pk,sf)
    Ny = Nx


    //---------------------------------------------------------------
    // Construction du maillage Pk pour pk = EF_Pk
    //---------------------------------------------------------------
    [Pk_Nv, Pk_Ne, Pk_CoorV, Pk_NuV, Pk_LogP, Pk_LogE, Pk_NuO] ...
                     = Struct_Pk_Mesh(Nx, Ny, Lx, Ly,EF_Pk)
    Pk_Nse = 3 ; 
    Pk_Nve = (EF_Pk+1)*(EF_Pk +2)/2.0

    // nombre  global de degré de liberté 
    Nv    = Pk_Nv   ;
    CoorV = Pk_CoorV ;
    LogP  = Pk_LogP ;

    // nombre local (sur chaque éléments) de degré de liberté 
    Nse = Pk_Nse;
    Nve = Pk_Nve;
    Ne  = Pk_Ne
    NuV = Pk_NuV

    //    disp(["rf, EF_Pk and Nv"])  
    //     disp([rf, EF_Pk, Nv])  

    // Initialisation de la matrice globale
    MatA(1:Nv,1:Nv)   = 0.0;

    // Initialisation du second membre gobal
    VecB(1:Nv)        = 0.0;

    //=========================================================
    // Assemblage 
    //=========================================================

    // boucle sur les éléments
    // ----------------------------------------------------
    for e =1:Ne

        Xi = CoorV(:, NuV(1,e));
        Xj = CoorV(:, NuV(2,e));
        Xk = CoorV(:, NuV(3,e));

        DetE = abs(Determinant(Xj-Xi,Xk-Xi)); 

        // Boucle sur les points de Gauss
        // ----------------------------------------------------------
        for gp = 1:Ng

            // Point physique associé au point de Gauss
            // ------------------------------------------
            Xg  = Lambda(1,gp)*Xi + Lambda(2,gp)*Xj + Lambda(3,gp)*Xk ;
            wg  = W(gp);

            [betx, bety , cx, cy, alpha] = Composite_Mat(Xg,choix)
            // betx=bety=dh
            //betx=sqrt(2.0/Ne)
            //bety=sqrt(2.0/Ne)
            
            // betx=bety=dh²
            //betx=(sqrt(2.0/Ne))*(sqrt(2.0/Ne))
            //bety=(sqrt(2.0/Ne))*(sqrt(2.0/Ne))
            
            // betx=bety=dh³
            //betx=(sqrt(2.0/Ne))*(sqrt(2.0/Ne))*(sqrt(2.0/Ne))
            //bety=(sqrt(2.0/Ne))*(sqrt(2.0/Ne))*(sqrt(2.0/Ne))

            GradP1             = MGrad_Lambda_of_P1(Xi, Xj, Xk);


            // Fonctions de base Pk
            // *****************************************
            Phi      = Phi_Pk(Lambda(1,gp), Lambda(2,gp), EF_Pk)

            // Gradient des fonctions de base Pk
            // *********************************************
            GradPhi  = GradPhi_Pk(Lambda(1,gp), Lambda(2,gp), ...
            GradP1(1,:),    GradP1(2,:),EF_Pk)

            // -----------------------------------------------------
            // Second membre local,  Matrice Locale
            //           et assemblage
            // -----------------------
            // Boucle sur la numérotation locale (lignes)
            for k=1:Nve
                // numérotation globale (is) associée à l'indice local (k)
                is          = NuV(k,e)      ;
                Phi_is      = Phi(k)        ; 
                GradPhi_is  = GradPhi(k, :) ;
                Be_k        = f(Xg, gamx, gamy, betx, bety, cx, cy, alpha)*Phi_is ;

                VecB(is)    = VecB(is) + wg*DetE*Be_k ;

                // Boucle sur la numérotation locale (colonnes)
                Beta = [betx, 0; 0, bety]; 
                C = [cx, cy]; 
                for kp = 1:Nve
                    // numérotation globale (js) à l'indice local (kp)
                    js         = NuV(kp,e)     ;
                    Phi_js     = Phi(kp)       ; 
                    GradPhi_js = GradPhi(kp,:) ;
                    Ae_k_kp    = (GradPhi_is * Beta) * GradPhi_js' + alpha * Phi_is * Phi_js + -Phi_is * (C * GradPhi_js') ;

                    MatA(is,js)= MatA(is,js) + wg*DetE*Ae_k_kp ;
                end // sur (kp)
            end  // sur (k) 
        end  // sur (g) 
    end  // sur (e)      
    
    // =====================================================
    // Conditions aux limites (peut être optimisé)
    //    Dirichlet homogène sur le bord gauche
    //    Neumann sur tous les autres bords
    // -------------------------------------------
    for is=1:Nv
        Xg = CoorV(:, is);

        if LogP(is) == -1 then
            MatA(is, :)    = 0.0
            MatA(is, is)   = 1.0
            VecB(is)       = Exact(Xg(1), Xg(2), gamx, gamy)
        end 
        if LogP(is) == -2 then
            MatA(is, :)  = 0.0
            MatA(is, is) = 1.0
            VecB(is)     = Exact(Xg(1), Xg(2), gamx, gamy)
        end 
        if LogP(is) < 0 then
            MatA(is, :)    = 0.0
            MatA(is, is)   = 1.0
            VecB(is)       = Exact(Xg(1), Xg(2), gamx, gamy)
        end     
    end


    // Version de résolution du système en 
    // utilisant une matrice creuse et umfpack
    // ---------------------------------------
    P1_Sparse_Mat = sparse(MatA); 
    VecSol        = umfpack(P1_Sparse_Mat,"\",VecB)
    //VecSol = MatA\VecB;

    //
    // Visualisation
    // -------------------------------------------
    Pk_x = CoorV(1,:)';
    Pk_y = CoorV(2,:)';


    je = 0
    for ie=1:Ne
        for jl = 1:int(size(Pk_NuO,2) )
            i1 = Pk_NuO(1,jl) ; 
            i2 = Pk_NuO(2,jl) ; 
            i3 = Pk_NuO(3,jl)  ;
            je = je +1
            Pk_P1_T(1:5,je) =[je,NuV(i1,ie),NuV(i2,ie), NuV(i3,ie),je] 
        end
    end


    Pk_exact        = Exact(Pk_x, Pk_y, gamx, gamy );
    gcf().color_map = jetcolormap(64);
    fec(Pk_x,Pk_y,Pk_P1_T',VecSol,mesh=%t)

    // boucle sur les éléments
    // ----------------------------------------------------
    // Calcul de l'erreur
    // ------------------------
   
    H(rf)   = sqrt(Lx*Ly/Ne)
    NvC(rf) = Nv; 
    VNx(rf) = Nx;

    ErrL1(rf)   = 0;   ErrL2(rf)   = 0; ErrLinf(rf) = 0.0 ;
    // boucle sur les éléments
    for e = 1:Ne
        // Boucle sur les points de Gauss
        Xi = CoorV(:, NuV(1,e));
        Xj = CoorV(:, NuV(2,e));
        Xk = CoorV(:, NuV(3,e));

        DetE          = abs(Determinant(Xi-Xj,Xk-Xj))      ;
        // Boucle sur les points de Gauss
        // ----------------------------------------------------------
        for gp = 1:Ng
            // Fonctions de base Pk
            // *****************************************
            Phi      = Phi_Pk(Lambda(1,gp), Lambda(2,gp), EF_Pk)
            wg   = W(gp);

            // Point physique associé au point de Gauss
            // ------------------------------------------        
            Xg(1:2) = 0 ;
            for k=1:Nse
                is      = NuV(k,e)     ;     
                Xg      = Xg + CoorV(:, is)*Lambda(k,gp); 
            end 

            // Solution approchée au point de Gauss
            // ------------------------------------------        
            Ug = 0
            for k=1:Nve
                is   = NuV(k,e)                 ;       
                Ug   = Ug +   VecSol(is)*Phi(k);
            end 

            // Solution exacte au point de Gauss
            // -----------------------------------------
            Uex  = Exact(Xg(1), Xg(2), gamx, gamy )    ;

            // Les Erreurs L1, L2 et Linf
            // ----------------------------------------------
            ErrL1(rf)   = ErrL1(rf) + wg*DetE*abs(Ug -Uex)  ; 
            ErrL2(rf)   = ErrL2(rf) + wg*DetE*(Ug -Uex)^2  ;
            ErrLinf(rf) = max( ErrLinf(rf), abs(Ug -Uex) ) ; 
        end     
    end
    ErrL2(rf) = sqrt(ErrL2(rf))

end

Dh=H(1:Nconv-1)./H(2:Nconv)
Ordre_de_Conv(1,:) = log(ErrL1(1:Nconv-1)./ErrL1(2:Nconv))./log(Dh)
Ordre_de_Conv(2,:) = log(ErrL2(1:Nconv-1)./ErrL2(2:Nconv))./log(Dh)
Ordre_de_Conv(3,:) = log(ErrLinf(1:Nconv-1)./ErrLinf(2:Nconv))./log(Dh)
disp("Elements Finis Pk avec k = ",EF_Pk)
disp("---------------------------------------------------------")
disp("choix = ",choix)
disp("---------------------------------------------------------")
disp("******* Ordre de convergence estimé *********")

disp([... 
"Nv",  "Norme L1","Norme L2","Norme Linf", "Erreur"] )
disp([... 
      NvC(2:Nconv)';...
      Ordre_de_Conv;ErrL2(2:Nconv)']')
disp("---------------------------------------------------------")

