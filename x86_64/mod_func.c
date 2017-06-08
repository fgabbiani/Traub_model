#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _traub_Ca_reg(void);
extern void _traub_Cashell_reg(void);
extern void _traub_Ka_reg(void);
extern void _traub_Kahp_reg(void);
extern void _traub_Kc_reg(void);
extern void _traub_Kdr_reg(void);
extern void _traub_Na_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," traub_Ca.mod");
    fprintf(stderr," traub_Cashell.mod");
    fprintf(stderr," traub_Ka.mod");
    fprintf(stderr," traub_Kahp.mod");
    fprintf(stderr," traub_Kc.mod");
    fprintf(stderr," traub_Kdr.mod");
    fprintf(stderr," traub_Na.mod");
    fprintf(stderr, "\n");
  }
  _traub_Ca_reg();
  _traub_Cashell_reg();
  _traub_Ka_reg();
  _traub_Kahp_reg();
  _traub_Kc_reg();
  _traub_Kdr_reg();
  _traub_Na_reg();
}
