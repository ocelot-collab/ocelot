/*
   Filename: eq_ens.h

   ENS server access class

   Created by: Arthur Aghababyan DESY-MVP

*/

#ifndef eq_ens_h
#define eq_ens_h

#include "eq_stl.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netdb.h>
#include "mtp.h"
#include "ens.h"
#include "eq_rpc.h"
#include "eq_adr.h"
#include "eq_data.h"

#ifdef __APPLE__
#include <sys/param.h>
#endif

#define DISCONNECT_ENS	   120


class eq_ens    {

      private:

        SEM_T                    sem;

        char                     host [MAXHOSTNAMELEN];

        char                     *list;
        char                     *evp;
        int                      reconnect;
        int                      num;

        time_t                   last_call;

        std::vector<std::string *> tine_context;

        CLIENT                   *cp;

        void    add_context      (EqAdr *);
        int     check_context    (EqAdr *);

        int     ens_num          (char *);
        int     ens_connect      (void);
        void    ens_disconnect   (void);
        void    ens_copy         (ENSinfo *, ENSinfo *);

        void    ens_error        (ENSinfo *);

        int     rpc_get          (ens_param *, ens_results *, int *);
        int     rpc_names        (ens_name_param *, ens_name_results *, int *);

        int     rpc_get_free     (ens_results *);
        int     rpc_names_free   (ens_name_results *);

        int     rpc_hard_error   (int);

        int     tine_type        (int, int);
        int     tine_get         (EqAdr *, ENSinfo *);
        int     tine_names       (EqAdr *, EqData *);

      public:
                eq_ens           (void);
                ~eq_ens          (void);

        int     ens_resolve      (EqAdr *, ENSinfo *);
        int     ens_names        (EqAdr *, EqData *);

        int     disconnect       (void);
};

#endif
