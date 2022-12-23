Uslip = - St
         *
         (I-St*Gr.T())
         &
        (
          fvc::ddt(U) +
          (U & fvc::grad(U)) +
          g
        );

phislip = fvc::flux(Uslip);
