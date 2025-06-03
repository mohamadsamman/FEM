clear 

// ///////////////////////////////////////////////////////////////////////////
//  Programme de calcul numérique sur un élément défini par x1, x2
//  l'intégration numérique est défini par 
//   Ng,  Lambda(1:Ng, 1) == Lambda_1  et Omega(1:Ng) 
//   on en déduit 
//    Lambda(1:Ng, 2)  = 1 -  Lambda(1:Ng, 1) == Lambda_2
// ----------------------------------------------------------------
//  les paramètres de la fonction numérique sont:
//       x1, x2
//       Ng, omega(1:Ng), Lambda1(1:Ng), Lambda2(1:Ng)
//       une fonction MyF(Lambda1, Lambda2)
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
function out = MyF(x)
    out = x
endfunction

function out = MyA11(L1,L2, dl1, dl2)
    out = dl1*dl1 + L1*L1 
endfunction
function out = MyA12(L1,L2, dl1, dl2)
    out = dl1*dl2 + L1*L2 
endfunction
function out = MyA21(L1,L2, dl1, dl2)
    out = dl2*dl1 + L2*L1
endfunction
function out = MyA22(L1,L2, dl1, dl2)
    out = dl2*dl2 + L2*L2 
endfunction

function out = MyF1(L1,L2,x)
    out = L1*MyF(x)
endfunction
function out = MyF2(L1,L2,x)
    out = L2*MyF(x)
endfunction

// --------------------------------------------------------------
// Intégration numérique d'une fonction 
// --------------------------------------------------------------
function out = IntegElmtA(x1,x2, Ng, Omega, Lambda1, Lambda2, dl1,dl2, MyFm)
    // Calcul de la somme qui définie l'intégration numérique.
    // -------------------------------------------------------
    out = 0 
    for g=1:Ng
        out = out + Omega(g)*MyFm(Lambda1(g), Lambda2(g), dl1, dl2 )
    end
    out = out*abs(x1 -x2)
endfunction
function out = IntegElmtF(x1,x2, Ng, Omega, Lambda1, Lambda2, MyFm)
    // Calcul de la somme qui définie l'intégration numérique.
    // -------------------------------------------------------
    out = 0 
    for g=1:Ng
        xg  = x1*Lambda1(g) + x2*Lambda2(g) 
        out = out + Omega(g)*MyFm(Lambda1(g), Lambda2(g),xg )
    end
    out = out*abs(x1 -x2)
endfunction


//--------------------------------------
// fonction second membre de problème.
//--------------------------------------
/* function  out = Rhs_Function(icas, x)
    select icas 
    case 1
        out =  x      //(%pi*a)*(%pi*a)*sin(%pi*a*x/Lx )/(Lx*Lx)
    case 2
        out =  x - 1  //(%pi*a)*(%pi*a)*sin(%pi*a*x/Lx )/(Lx*Lx)
    else
        out =  x      //(%pi*a)*(%pi*a)*sin(%pi*a*x/Lx )/(Lx*Lx)
    end
    
endfunction
*/

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
    
// fonctions de base aux points de Gauss
// ====================================================
x1e = 2 ; x2e = 3  ; Nve = 3

Grad_Lambda = MGrad_Lambda_of_P1(x1e, x2e)
dl1dx = Grad_Lambda(1);
dl2dx = Grad_Lambda(2);

AeP1(1,1) = IntegElmtA(x1e,x2e, Ng, Omega, Lambda(1,:), Lambda(2,:), ...
                           dl1dx, dl2dx, MyA11)
AeP1(1,2) = IntegElmtA(x1e,x2e, Ng, Omega, Lambda(1,:), Lambda(2,:), ...
                           dl1dx, dl2dx, MyA12)
AeP1(2,1) = IntegElmtA(x1e,x2e, Ng, Omega, Lambda(1,:), Lambda(2,:), ...
                           dl1dx, dl2dx, MyA21)
AeP1(2,2) = IntegElmtA(x1e,x2e, Ng, Omega, Lambda(1,:), Lambda(2,:), ...
                           dl1dx, dl2dx, MyA22)
                           
BeP1(1)   = IntegElmtF(x1e,x2e, Ng, Omega, Lambda(1,:), Lambda(2,:), MyF1)
BeP1(2)   = IntegElmtF(x1e,x2e, Ng, Omega, Lambda(1,:), Lambda(2,:), MyF2)

// Solution calculée analytiquement :: Exercice 2.1.1
Aex=[4/3, -5/6 ; -5/6, 4/3]
Bex=[7/6;8/6]

disp(norm(AeP1-Aex))
disp(norm(BeP1-Bex))
/*

// fonctions de base et gradients aux points de Gauss
// ====================================================
select Nve 
case 2
        Pk_Bf_at_Xg(1,1:Ng) =  Lambda(1,1:Ng)   
        Pk_Bf_at_Xg(2,1:Ng) =  1.0 - Lambda(1, 1:Ng); 

       Grad_Lambda = MGrad_Lambda_of_P1(x1e, x2e) 
           
       Grad_Pk_Bf_at_Xg(1,1:Ng) = Grad_Lambda(1)
       Grad_Pk_Bf_at_Xg(2,1:Ng) = Grad_Lambda(2)
       
       Aex=[4/3, -5/6 ; -5/6, 4/3]
       Bex=[7/6;8/6]

case 3
        Pk_Bf_at_Xg(1,1:Ng) =  Lambda(1,1:Ng).*( 2.0*Lambda(1,1:Ng) - 1)   
        Pk_Bf_at_Xg(2,1:Ng) =  Lambda(2,1:Ng).*( 2.0*Lambda(2,1:Ng) - 1)
        Pk_Bf_at_Xg(3,1:Ng) =  4.0*Lambda(1,1:Ng).*Lambda(2,1:Ng)
        
       Grad_Lambda = MGrad_Lambda_of_P1(x1e, x2e) 
           
       Grad_Pk_Bf_at_Xg(1,1:Ng) = Grad_Lambda(1)*(4.0*Lambda(1,1:Ng) - 1)
       Grad_Pk_Bf_at_Xg(2,1:Ng) = Grad_Lambda(2)*(4.0*Lambda(2,1:Ng) - 1) 
       Grad_Pk_Bf_at_Xg(3,1:Ng) = -4*Grad_Lambda(2)*(2.0*Lambda(2,1:Ng) - 1)       
       Aex=[7/3, 1/3, -8/3; 1/3, 7/3, -8/3; -8/3, -8/3, 16/3]
      Bex=[1/6;2/6;1]

end

Icas = 2
select Icas 
case 1
Coef_Mu = 1.0
Coef_Cv = 1.0

case 2
Coef_Mu = 1.0
Coef_Cv = 0.0    
end

Mat_Ae(1:Nve,1:Nve) = 0
Vec_Be(1:Nve)       = 0

aire         = abs(x2e-x1e);

// boucles sur les points de Gauss
// ===============================
for g = 1:Ng

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

        fg        = Rhs_Function(Icas, xg)

        bl        = aire*fg*Vl 

        Vec_Be(l) =  Vec_Be(l)  + Omega(g)*bl 
        

        // seconde boucle (inbriquée) sur les ddl locaux
        // =============================================
        // Assemblage local de la matrice 
        for m = 1:Nve
            Vm       =      Pk_Bf_at_Xg(m,g)
            GVm      = Grad_Pk_Bf_at_Xg(m,g)
            
            Alm       = aire * ( Coef_Mu*GVl*GVm + Coef_Cv*Vl*Vm )

            Mat_Ae(l,m) =  Mat_Ae(l,m) + Omega(g)*Alm 
        end
    end
end
// Conditions aux limites de Dirichlet homogène
Ibord = 2
select Ibord
case 1
 // Conditions aux limites de Dirichlet homogène en l= 1 et l=2
 Mat_Ae(1,:) = 0 ; Mat_Ae(1,1) = 1 ; Vec_Be(1) = 0
 Mat_Ae(2,:) = 0 ; Mat_Ae(2,2) = 1 ; Vec_Be(2) = 0  
 Exact=[0;0;3/16]
case 2
 // Conditions de Dirichlet homogène en  l=2
 Mat_Ae(2,:) = 0 ; Mat_Ae(2,2) = 1 ; Vec_Be(2) = 0
 Exact=[2/3;0;25/48]
end

Resu = Mat_Ae\Vec_Be
// disp(Resu)
disp(norm(Resu-Exact))

*/
