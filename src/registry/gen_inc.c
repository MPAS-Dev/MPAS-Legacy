#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dictionary.h"
#include "registry_types.h"
#include "gen_inc.h"
#include "fortprintf.h"

int is_derived_dim(char * d)
{
   if (strchr(d, (int)'+')) return 1;
   if (strchr(d, (int)'-')) return 1;

   return 0;
}

void split_derived_dim_string(char * dim, char ** p1, char ** p2)
{
   char * cp, * cm, * c;
   int n;

   cp = strchr(dim, (int)'+');
   cm = strchr(dim, (int)'-');
   if (!cp) 
      c = cm;
   else if (!cm) 
      c = cp;
   else if (cm < cp) 
      c = cm;
   else 
      c = cp;

   n = c - dim;
   *p1 = (char *)malloc(n*sizeof(char));
   snprintf(*p1, n, "%s", dim+1);

   *p2 = (char *)malloc((strlen(dim)-n+1)*sizeof(char));
   sprintf(*p2, "%s", dim+n);
}

void gen_namelists(struct namelist * nls)
{
   struct namelist * nls_ptr;
   struct dtable * dictionary;
   int done;
   char nlrecord[1024];
   FILE * fd;

   /*
    *  Generate config_type.inc
    */
   fd = fopen("config_defs.inc", "w");

   nls_ptr = nls;
   while (nls_ptr) {
      if (nls_ptr->vtype == INTEGER)   fortprintf(fd, "   integer :: %s\n",nls_ptr->name);
      if (nls_ptr->vtype == REAL)      fortprintf(fd, "   real (KIND=RKIND) :: %s\n",nls_ptr->name);
      if (nls_ptr->vtype == LOGICAL)   fortprintf(fd, "   logical :: %s\n",nls_ptr->name);
      if (nls_ptr->vtype == CHARACTER) fortprintf(fd, "   character (len=32) :: %s\n",nls_ptr->name);

      nls_ptr = nls_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate namelist_defs.inc
    */
   fd = fopen("config_namelist_defs.inc", "w");
   dict_alloc(&dictionary);

   done = 0;
  
   while (!done) {
      nls_ptr = nls;
      while (nls_ptr && dict_search(dictionary, nls_ptr->record))
         nls_ptr = nls_ptr->next;

      if (nls_ptr) {
         dict_insert(dictionary, nls_ptr->record);
         strncpy(nlrecord, nls_ptr->record, 1024);
         fortprintf(fd, "      namelist /%s/ %s", nls_ptr->record, nls_ptr->name);
         nls_ptr = nls_ptr->next;
         while(nls_ptr) {
            if (strncmp(nls_ptr->record, nlrecord, 1024) == 0)
               fortprintf(fd, ", &\n                    %s", nls_ptr->name);
            nls_ptr = nls_ptr->next;
         }
         fortprintf(fd, "\n");
      }
      else
         done = 1;
   }
   

   dict_free(&dictionary);
   fclose(fd);


   /*
    *  Generate namelist_reads.inc
    */
   fd = fopen("config_set_defaults.inc", "w");
   nls_ptr = nls;
   while (nls_ptr) {
      if (nls_ptr->vtype == INTEGER) fortprintf(fd, "      %s = %i\n", nls_ptr->name, nls_ptr->defval.ival);
      if (nls_ptr->vtype == REAL)    fortprintf(fd, "      %s = %f\n", nls_ptr->name, nls_ptr->defval.rval);
      if (nls_ptr->vtype == LOGICAL) {
         if (nls_ptr->defval.lval == 0) 
            fortprintf(fd, "      %s = .false.\n", nls_ptr->name);
         else
            fortprintf(fd, "      %s = .true.\n", nls_ptr->name);
      }
      if (nls_ptr->vtype == CHARACTER)
         fortprintf(fd, "      %s = \"%s\"\n", nls_ptr->name, nls_ptr->defval.cval);
      nls_ptr = nls_ptr->next;
   }
   fortprintf(fd, "\n");
   fclose(fd);


   fd = fopen("config_namelist_reads.inc", "w");
   dict_alloc(&dictionary);
   nls_ptr = nls;
   while (nls_ptr) {
      if (!dict_search(dictionary, nls_ptr->record)) {
         fortprintf(fd, "         read(funit,%s)\n", nls_ptr->record);
         dict_insert(dictionary, nls_ptr->record);
      }
      nls_ptr = nls_ptr->next;
   }
   fortprintf(fd, "\n");

   dict_free(&dictionary);
   fclose(fd);


   fd = fopen("config_bcast_namelist.inc", "w");
   nls_ptr = nls;
   while (nls_ptr) {
      if (nls_ptr->vtype == INTEGER)   fortprintf(fd, "      call dmpar_bcast_int(dminfo, %s)\n", nls_ptr->name);
      if (nls_ptr->vtype == REAL)      fortprintf(fd, "      call dmpar_bcast_real(dminfo, %s)\n", nls_ptr->name);
      if (nls_ptr->vtype == LOGICAL)   fortprintf(fd, "      call dmpar_bcast_logical(dminfo, %s)\n", nls_ptr->name);
      if (nls_ptr->vtype == CHARACTER) fortprintf(fd, "      call dmpar_bcast_char(dminfo, %s)\n", nls_ptr->name);
      nls_ptr = nls_ptr->next;
   }
   fortprintf(fd, "\n");
   fclose(fd);
}


void gen_field_defs(struct variable * vars, struct dimension * dims)
{
   struct variable * var_ptr;
   struct variable * var_ptr2;
   struct dimension * dim_ptr;
   struct dimension_list * dimlist_ptr;
   FILE * fd;
   char super_array[1024];
   char array_class[1024];
   int i;
   int class_start, class_end;
   int vtype;

   /*
    * Generate indices for super arrays
    */
   fd = fopen("super_array_indices.inc", "w");
   var_ptr = vars;
   memcpy(super_array, var_ptr->super_array, 1024);
   i = 1;
   while (var_ptr) {
      if (strncmp(super_array, var_ptr->super_array, 1024) != 0) {
         memcpy(super_array, var_ptr->super_array, 1024);
         i = 1;
       }
      if (strncmp(var_ptr->array_class, "-", 1024) != 0) fortprintf(fd, "      integer :: index_%s = %i\n", var_ptr->name_in_code, i++);
      var_ptr = var_ptr->next;
   }
   var_ptr = vars;
   memcpy(super_array, var_ptr->super_array, 1024);
   memcpy(array_class, var_ptr->array_class, 1024);
   class_start = 1;
   class_end = 1;
   i = 1;
   while (var_ptr) {
      if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
         if (strncmp(super_array, var_ptr->super_array, 1024) != 0) {
            if (strncmp(super_array, "-", 1024) != 0) fortprintf(fd, "      integer :: %s_end = %i\n", array_class, class_end);
            if (strncmp(super_array, "-", 1024) != 0) fortprintf(fd, "      integer :: num_%s = %i\n", super_array, i);
            class_start = 1;
            class_end = 1;
            i = 1;
            memcpy(super_array, var_ptr->super_array, 1024);
            memcpy(array_class, var_ptr->array_class, 1024);
            fortprintf(fd, "      integer :: %s_start = %i\n", array_class, class_start);
         }
         else if (strncmp(array_class, var_ptr->array_class, 1024) != 0) {
            fortprintf(fd, "      integer :: %s_end = %i\n", array_class, class_end);
            class_start = class_end+1;
            class_end = class_start;
            memcpy(array_class, var_ptr->array_class, 1024);
            fortprintf(fd, "      integer :: %s_start = %i\n", array_class, class_start);
            i++;
         }
         else {
            class_end++;
            i++;
         }
      }
      var_ptr = var_ptr->next;
   }
   if (strncmp(super_array, "-", 1024) != 0) fortprintf(fd, "      integer :: %s_end = %i\n", array_class, class_end);
   if (strncmp(super_array, "-", 1024) != 0) fortprintf(fd, "      integer :: num_%s = %i\n", super_array, i);
   fclose(fd);


   /*
    *  Generate declarations of dimensions
    */
   fd = fopen("field_dimensions.inc", "w");
   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: %s\n", dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: %s\n", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: %sSolve\n", dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: %sSolve\n", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate dummy dimension argument list
    */
   fd = fopen("dim_dummy_args.inc", "w");
   dim_ptr = dims;
   if (dim_ptr && dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
      fortprintf(fd, "                            %s", dim_ptr->name_in_code);
      dim_ptr = dim_ptr->next;
   }
   else if (dim_ptr && dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
      fortprintf(fd, "                            %s", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, ", %s", dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, ", %s", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, " &\n");

   fclose(fd);


   /*
    *  Generate dummy dimension argument declaration list
    */
   fd = fopen("dim_dummy_decls.inc", "w");
   dim_ptr = dims;
   if (dim_ptr && dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
      fortprintf(fd, "      integer, intent(in) :: %s", dim_ptr->name_in_code);
      dim_ptr = dim_ptr->next;
   }
   else if (dim_ptr && dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
      fortprintf(fd, "      integer, intent(in) :: %s", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, ", %s", dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, ", %s", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   fclose(fd);


   /*
    *  Generate declarations of dimensions
    */
   fd = fopen("dim_decls.inc", "w");
   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: %s\n", dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: %s\n", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate calls to read dimensions from input file
    */
   fd = fopen("read_dims.inc", "w");
   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      call io_input_get_dimension(input_obj, \'%s\', %s)\n", dim_ptr->name_in_file, dim_ptr->name_in_code);
      else if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      call io_input_get_dimension(input_obj, \'%s\', %s)\n", dim_ptr->name_in_file, dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate declarations of time-invariant fields
    */
   fd = fopen("time_invariant_fields.inc", "w");
   var_ptr = vars;
   while (var_ptr) {
      if (var_ptr->timedim == 0) {
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            memcpy(super_array, var_ptr->super_array, 1024);
            memcpy(array_class, var_ptr->array_class, 1024);
            vtype = var_ptr->vtype;
            while (var_ptr && strncmp(super_array, var_ptr->super_array, 1024) == 0) {
               var_ptr2 = var_ptr;
               var_ptr = var_ptr->next;
            }
            if (vtype == INTEGER)  fortprintf(fd, "      type (field%idInteger), pointer :: %s\n", (var_ptr2->ndims)+1, var_ptr2->super_array);
            if (vtype == REAL)     fortprintf(fd, "      type (field%idReal), pointer :: %s\n", (var_ptr2->ndims)+1, var_ptr2->super_array);
         }
         else {
            if (var_ptr->vtype == INTEGER)  fortprintf(fd, "      type (field%idInteger), pointer :: %s\n", var_ptr->ndims, var_ptr->name_in_code);
            if (var_ptr->vtype == REAL)     fortprintf(fd, "      type (field%idReal), pointer :: %s\n", var_ptr->ndims, var_ptr->name_in_code);
            var_ptr = var_ptr->next;
         }
      }
      else
         var_ptr = var_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate declarations of time-invariant fields
    */
   fd = fopen("time_varying_fields.inc", "w");
   var_ptr = vars;
   while (var_ptr) {
      if (var_ptr->timedim == 1) {
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            memcpy(super_array, var_ptr->super_array, 1024);
            memcpy(array_class, var_ptr->array_class, 1024);
            vtype = var_ptr->vtype;
            while (var_ptr && strncmp(super_array, var_ptr->super_array, 1024) == 0) {
               var_ptr2 = var_ptr;
               var_ptr = var_ptr->next;
            }
            if (vtype == INTEGER)  fortprintf(fd, "      type (field%idInteger), pointer :: %s\n", (var_ptr2->ndims)+1, var_ptr2->super_array);
            if (vtype == REAL)     fortprintf(fd, "      type (field%idReal), pointer :: %s\n", (var_ptr2->ndims)+1, var_ptr2->super_array);
         }
         else {
            if (var_ptr->vtype == INTEGER)  fortprintf(fd, "      type (field%idInteger), pointer :: %s\n", var_ptr->ndims, var_ptr->name_in_code);
            if (var_ptr->vtype == REAL)     fortprintf(fd, "      type (field%idReal), pointer :: %s\n", var_ptr->ndims, var_ptr->name_in_code);
            var_ptr = var_ptr->next;
         }
      }
      else
         var_ptr = var_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate grid metadata allocations
    */
   fd = fopen("grid_meta_allocs.inc", "w");

   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      g %% %s = %s\n", dim_ptr->name_in_code, dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      g %% %s = %s\n", dim_ptr->name_in_file, dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   var_ptr = vars;
   while (var_ptr) {
      if (var_ptr->timedim == 0) {
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            memcpy(super_array, var_ptr->super_array, 1024);
            memcpy(array_class, var_ptr->array_class, 1024);
            vtype = var_ptr->vtype;
            i = 0;
            while (var_ptr && strncmp(super_array, var_ptr->super_array, 1024) == 0) {
               i++;
               var_ptr2 = var_ptr;
               var_ptr = var_ptr->next;
            }
            fortprintf(fd, "      allocate(g %% %s)\n", var_ptr2->super_array);
            fortprintf(fd, "      allocate(g %% %s %% ioinfo)\n", var_ptr2->super_array);
            fortprintf(fd, "      allocate(g %% %s %% array(%i, ", var_ptr2->super_array, i);
            dimlist_ptr = var_ptr2->dimlist;
            if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
               fortprintf(fd, "%s + 1", dimlist_ptr->dim->name_in_code);
            else
               if (dimlist_ptr->dim->namelist_defined) fortprintf(fd, "%s", dimlist_ptr->dim->name_in_file);
               else fortprintf(fd, "%s", dimlist_ptr->dim->name_in_code);
            dimlist_ptr = dimlist_ptr->next;
            while (dimlist_ptr) {
               if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
                  fortprintf(fd, ", %s + 1", dimlist_ptr->dim->name_in_code);
               else
                  if (dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_file);
                  else fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_code);
               dimlist_ptr = dimlist_ptr->next;
            }
            fortprintf(fd, "))\n");
            fortprintf(fd, "      g %% %s %% array = 0\n", var_ptr2->super_array ); /* initialize field to zero */

            if (var_ptr2->iostreams & INPUT0) 
               fortprintf(fd, "      g %% %s %% ioinfo %% input = .true.\n", var_ptr2->super_array);
            else
               fortprintf(fd, "      g %% %s %% ioinfo %% input = .false.\n", var_ptr2->super_array);

            if (var_ptr2->iostreams & RESTART0) 
               fortprintf(fd, "      g %% %s %% ioinfo %% restart = .true.\n", var_ptr2->super_array);
            else
               fortprintf(fd, "      g %% %s %% ioinfo %% restart = .false.\n", var_ptr2->super_array);

            if (var_ptr2->iostreams & OUTPUT0) 
               fortprintf(fd, "      g %% %s %% ioinfo %% output = .true.\n", var_ptr2->super_array);
            else
               fortprintf(fd, "      g %% %s %% ioinfo %% output = .false.\n", var_ptr2->super_array);
            fortprintf(fd, "\n");
         }
         else {
            fortprintf(fd, "      allocate(g %% %s)\n", var_ptr->name_in_code);
            fortprintf(fd, "      allocate(g %% %s %% ioinfo)\n", var_ptr->name_in_code);
            if (var_ptr->ndims > 0) {
               fortprintf(fd, "      allocate(g %% %s %% array(", var_ptr->name_in_code);
               dimlist_ptr = var_ptr->dimlist;
               if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
                  fortprintf(fd, "%s + 1", dimlist_ptr->dim->name_in_code);
               else
                  if (dimlist_ptr->dim->namelist_defined) fortprintf(fd, "%s", dimlist_ptr->dim->name_in_file);
                  else fortprintf(fd, "%s", dimlist_ptr->dim->name_in_code);
               dimlist_ptr = dimlist_ptr->next;
               while (dimlist_ptr) {
                  if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                      !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                      !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
                     fortprintf(fd, ", %s + 1", dimlist_ptr->dim->name_in_code);
                  else
                     if (dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_file);
                     else fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_code);
                  dimlist_ptr = dimlist_ptr->next;
               }
               fortprintf(fd, "))\n");
               fortprintf(fd, "      g %% %s %% array = 0\n", var_ptr->name_in_code ); /* initialize field to zero */

            }
            if (var_ptr->iostreams & INPUT0) 
               fortprintf(fd, "      g %% %s %% ioinfo %% input = .true.\n", var_ptr->name_in_code);
            else
               fortprintf(fd, "      g %% %s %% ioinfo %% input = .false.\n", var_ptr->name_in_code);

            if (var_ptr->iostreams & RESTART0) 
               fortprintf(fd, "      g %% %s %% ioinfo %% restart = .true.\n", var_ptr->name_in_code);
            else
               fortprintf(fd, "      g %% %s %% ioinfo %% restart = .false.\n", var_ptr->name_in_code);

            if (var_ptr->iostreams & OUTPUT0) 
               fortprintf(fd, "      g %% %s %% ioinfo %% output = .true.\n", var_ptr->name_in_code);
            else
               fortprintf(fd, "      g %% %s %% ioinfo %% output = .false.\n", var_ptr->name_in_code);
            fortprintf(fd, "\n");

            var_ptr = var_ptr->next;
         }
      }
      else
         var_ptr = var_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate grid metadata deallocations
    */
   fd = fopen("grid_meta_deallocs.inc", "w");

   var_ptr = vars;
   while (var_ptr) {
      if (var_ptr->timedim == 0) {
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            memcpy(super_array, var_ptr->super_array, 1024);
            memcpy(array_class, var_ptr->array_class, 1024);
            vtype = var_ptr->vtype;
            i = 0;
            while (var_ptr && strncmp(super_array, var_ptr->super_array, 1024) == 0) {
               i++;
               var_ptr2 = var_ptr;
               var_ptr = var_ptr->next;
            }
            fortprintf(fd, "      deallocate(g %% %s %% array)\n", var_ptr2->super_array);
            fortprintf(fd, "      deallocate(g %% %s %% ioinfo)\n", var_ptr2->super_array);
            fortprintf(fd, "      deallocate(g %% %s)\n\n", var_ptr2->super_array);
         }
         else {
            if (var_ptr->ndims > 0) {
               fortprintf(fd, "      deallocate(g %% %s %% array)\n", var_ptr->name_in_code);
               fortprintf(fd, "      deallocate(g %% %s %% ioinfo)\n", var_ptr->name_in_code);
               fortprintf(fd, "      deallocate(g %% %s)\n\n", var_ptr->name_in_code);
            }
            else {
               fortprintf(fd, "      deallocate(g %% %s %% ioinfo)\n", var_ptr->name_in_code);
               fortprintf(fd, "      deallocate(g %% %s)\n\n", var_ptr->name_in_code);
            }
            var_ptr = var_ptr->next;
         }
      }
      else
         var_ptr = var_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate grid state allocations
    */
   fd = fopen("grid_state_allocs.inc", "w");

   var_ptr = vars;
   while (var_ptr) {
      if (var_ptr->timedim == 1 && var_ptr->ndims > 0) {
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            memcpy(super_array, var_ptr->super_array, 1024);
            memcpy(array_class, var_ptr->array_class, 1024);
            vtype = var_ptr->vtype;
            i = 0;
            while (var_ptr && strncmp(super_array, var_ptr->super_array, 1024) == 0) {
               i++;
               var_ptr2 = var_ptr;
               var_ptr = var_ptr->next;
            }
            fortprintf(fd, "      allocate(s %% %s)\n", var_ptr2->super_array);
            fortprintf(fd, "      allocate(s %% %s %% ioinfo)\n", var_ptr2->super_array);
            fortprintf(fd, "      allocate(s %% %s %% array(%i, ", var_ptr2->super_array, i);
            dimlist_ptr = var_ptr2->dimlist;
            if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
               fortprintf(fd, "b %% mesh %% %s + 1", dimlist_ptr->dim->name_in_code);
            else
               if (dimlist_ptr->dim->namelist_defined) fortprintf(fd, "b %% mesh %% %s", dimlist_ptr->dim->name_in_file);
               else fortprintf(fd, "b %% mesh %% %s", dimlist_ptr->dim->name_in_code);
            dimlist_ptr = dimlist_ptr->next;
            while (dimlist_ptr) {
               if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
                  fortprintf(fd, ", b %% mesh %% %s + 1", dimlist_ptr->dim->name_in_code);
               else
                  if (dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", b %% mesh %% %s", dimlist_ptr->dim->name_in_file);
                  else fortprintf(fd, ", b %% mesh %% %s", dimlist_ptr->dim->name_in_code);
               dimlist_ptr = dimlist_ptr->next;
            }
            fortprintf(fd, "))\n");
            fortprintf(fd, "      s %% %s %% block => b\n", var_ptr2->super_array);
   
            if (var_ptr2->iostreams & INPUT0) 
               fortprintf(fd, "      s %% %s %% ioinfo %% input = .true.\n", var_ptr2->super_array);
            else
               fortprintf(fd, "      s %% %s %% ioinfo %% input = .false.\n", var_ptr2->super_array);
   
            if (var_ptr2->iostreams & RESTART0) 
               fortprintf(fd, "      s %% %s %% ioinfo %% restart = .true.\n", var_ptr2->super_array);
            else
               fortprintf(fd, "      s %% %s %% ioinfo %% restart = .false.\n", var_ptr2->super_array);
   
            if (var_ptr2->iostreams & OUTPUT0) 
               fortprintf(fd, "      s %% %s %% ioinfo %% output = .true.\n", var_ptr2->super_array);
            else
               fortprintf(fd, "      s %% %s %% ioinfo %% output = .false.\n", var_ptr2->super_array);
            fortprintf(fd, "\n");
         }
         else {
            fortprintf(fd, "      allocate(s %% %s)\n", var_ptr->name_in_code);
            fortprintf(fd, "      allocate(s %% %s %% ioinfo)\n", var_ptr->name_in_code);
            fortprintf(fd, "      allocate(s %% %s %% array(", var_ptr->name_in_code);
            dimlist_ptr = var_ptr->dimlist;
            if (dimlist_ptr->dim->constant_value < 0) {
               if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
                  fortprintf(fd, "b %% mesh %% %s + 1", dimlist_ptr->dim->name_in_code);
               else
                  if (dimlist_ptr->dim->namelist_defined) fortprintf(fd, "b %% mesh %% %s", dimlist_ptr->dim->name_in_file);
                  else fortprintf(fd, "b %% mesh %% %s", dimlist_ptr->dim->name_in_code);
            }
            else {
               fortprintf(fd, "%i", dimlist_ptr->dim->constant_value);
            }
            dimlist_ptr = dimlist_ptr->next;
            while (dimlist_ptr) {
               if (dimlist_ptr->dim->constant_value < 0) {
                  if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                      !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                      !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
                     fortprintf(fd, ", b %% mesh %% %s + 1", dimlist_ptr->dim->name_in_code);
                  else
                     if (dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", b %% mesh %% %s", dimlist_ptr->dim->name_in_file);
                     else fortprintf(fd, ", b %% mesh %% %s", dimlist_ptr->dim->name_in_code);
               }
               else {
                  fortprintf(fd, ", %i", dimlist_ptr->dim->constant_value);
               }
               dimlist_ptr = dimlist_ptr->next;
            }
            fortprintf(fd, "))\n");
            fortprintf(fd, "      s %% %s %% block => b\n", var_ptr->name_in_code);
   
            if (var_ptr->iostreams & INPUT0) 
               fortprintf(fd, "      s %% %s %% ioinfo %% input = .true.\n", var_ptr->name_in_code);
            else
               fortprintf(fd, "      s %% %s %% ioinfo %% input = .false.\n", var_ptr->name_in_code);
   
            if (var_ptr->iostreams & RESTART0) 
               fortprintf(fd, "      s %% %s %% ioinfo %% restart = .true.\n", var_ptr->name_in_code);
            else
               fortprintf(fd, "      s %% %s %% ioinfo %% restart = .false.\n", var_ptr->name_in_code);
   
            if (var_ptr->iostreams & OUTPUT0) 
               fortprintf(fd, "      s %% %s %% ioinfo %% output = .true.\n", var_ptr->name_in_code);
            else
               fortprintf(fd, "      s %% %s %% ioinfo %% output = .false.\n", var_ptr->name_in_code);
            fortprintf(fd, "\n");
            var_ptr = var_ptr->next;
         }
      }
      else if (var_ptr->timedim == 1) {
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            memcpy(super_array, var_ptr->super_array, 1024);
            memcpy(array_class, var_ptr->array_class, 1024);
            vtype = var_ptr->vtype;
            i = 0;
            while (var_ptr && strncmp(super_array, var_ptr->super_array, 1024) == 0) {
               i++;
               var_ptr2 = var_ptr;
               var_ptr = var_ptr->next;
            }
            fortprintf(fd, "      allocate(s %% %s)\n", var_ptr2->super_array);
            fortprintf(fd, "      allocate(s %% %s %% ioinfo)\n", var_ptr2->super_array);
            fortprintf(fd, "      allocate(s %% %s %% array(%i)", var_ptr->name_in_code, i);
            fortprintf(fd, "      s %% %s %% block => b\n", var_ptr2->super_array);
   
            if (var_ptr2->iostreams & INPUT0) 
               fortprintf(fd, "      s %% %s %% ioinfo %% input = .true.\n", var_ptr2->super_array);
            else
               fortprintf(fd, "      s %% %s %% ioinfo %% input = .false.\n", var_ptr2->super_array);
   
            if (var_ptr2->iostreams & RESTART0) 
               fortprintf(fd, "      s %% %s %% ioinfo %% restart = .true.\n", var_ptr2->super_array);
            else
               fortprintf(fd, "      s %% %s %% ioinfo %% restart = .false.\n", var_ptr2->super_array);
   
            if (var_ptr2->iostreams & OUTPUT0) 
               fortprintf(fd, "      s %% %s %% ioinfo %% output = .true.\n", var_ptr2->super_array);
            else
               fortprintf(fd, "      s %% %s %% ioinfo %% output = .false.\n", var_ptr2->super_array);
            fortprintf(fd, "\n");
         }
         else {
            fortprintf(fd, "      allocate(s %% %s)\n", var_ptr->name_in_code);
            fortprintf(fd, "      allocate(s %% %s %% ioinfo)\n", var_ptr->name_in_code);
            fortprintf(fd, "      s %% %s %% block => b\n", var_ptr->name_in_code);
   
            if (var_ptr->iostreams & INPUT0) 
               fortprintf(fd, "      s %% %s %% ioinfo %% input = .true.\n", var_ptr->name_in_code);
            else
               fortprintf(fd, "      s %% %s %% ioinfo %% input = .false.\n", var_ptr->name_in_code);
   
            if (var_ptr->iostreams & RESTART0) 
               fortprintf(fd, "      s %% %s %% ioinfo %% restart = .true.\n", var_ptr->name_in_code);
            else
               fortprintf(fd, "      s %% %s %% ioinfo %% restart = .false.\n", var_ptr->name_in_code);
   
            if (var_ptr->iostreams & OUTPUT0) 
               fortprintf(fd, "      s %% %s %% ioinfo %% output = .true.\n", var_ptr->name_in_code);
            else
               fortprintf(fd, "      s %% %s %% ioinfo %% output = .false.\n", var_ptr->name_in_code);
            fortprintf(fd, "\n");
            var_ptr = var_ptr->next;
         }
      } 
      else
         var_ptr = var_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate grid state deallocations
    */
   fd = fopen("grid_state_deallocs.inc", "w");

   var_ptr = vars;
   while (var_ptr) {
      if (var_ptr->timedim == 1) {
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            memcpy(super_array, var_ptr->super_array, 1024);
            memcpy(array_class, var_ptr->array_class, 1024);
            vtype = var_ptr->vtype;
            i = 0;
            while (var_ptr && strncmp(super_array, var_ptr->super_array, 1024) == 0) {
               i++;
               var_ptr2 = var_ptr;
               var_ptr = var_ptr->next;
            }
            fortprintf(fd, "      deallocate(s %% %s %% array)\n", var_ptr2->super_array);
            fortprintf(fd, "      deallocate(s %% %s %% ioinfo)\n", var_ptr2->super_array);
            fortprintf(fd, "      deallocate(s %% %s)\n\n", var_ptr2->super_array);
         }
         else {
            if (var_ptr->ndims > 0) fortprintf(fd, "      deallocate(s %% %s %% array)\n", var_ptr->name_in_code);
            fortprintf(fd, "      deallocate(s %% %s %% ioinfo)\n", var_ptr->name_in_code);
            fortprintf(fd, "      deallocate(s %% %s)\n\n", var_ptr->name_in_code);
            var_ptr = var_ptr->next;
         }
      }
      else
         var_ptr = var_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate copies of state arrays
    */
   fd = fopen("copy_state.inc", "w");

   var_ptr = vars;
   while (var_ptr) {
      if (var_ptr->timedim == 1) {
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            memcpy(super_array, var_ptr->super_array, 1024);
            memcpy(array_class, var_ptr->array_class, 1024);
            vtype = var_ptr->vtype;
            i = 0;
            while (var_ptr && strncmp(super_array, var_ptr->super_array, 1024) == 0) {
               i++;
               var_ptr2 = var_ptr;
               var_ptr = var_ptr->next;
            }
            if (var_ptr2->ndims > 0) 
               fortprintf(fd, "      dest %% %s %% array = src %% %s %% array\n", var_ptr2->super_array, var_ptr2->super_array);
            else
               fortprintf(fd, "      dest %% %s %% scalar = src %% %s %% scalar\n", var_ptr2->super_array, var_ptr2->super_array);
         }
         else {
            if (var_ptr->ndims > 0) 
               fortprintf(fd, "      dest %% %s %% array = src %% %s %% array\n", var_ptr->name_in_code, var_ptr->name_in_code);
            else
               fortprintf(fd, "      dest %% %s %% scalar = src %% %s %% scalar\n", var_ptr->name_in_code, var_ptr->name_in_code);
            var_ptr = var_ptr->next;
         }
      }
      else
         var_ptr = var_ptr->next;
   }

   fclose(fd);
}


void gen_reads(struct variable * vars, struct dimension * dims)
{
   struct variable * var_ptr;
   struct dimension * dim_ptr;
   struct dimension_list * dimlist_ptr, * lastdim;
   struct dtable * dictionary;
   FILE * fd;
   char vtype[5];
   char fname[32];
   char * cp1, * cp2;
   int i, j;
   int ivtype;
   int has_vert_dim, vert_dim;


   /*
    *  Generate declarations of IDs belonging in io_input_object
    */
   fd = fopen("io_input_obj_decls.inc", "w");

   fortprintf(fd, "      integer :: rdDimIDTime\n");
   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: rdDimID%s\n", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   fortprintf(fd, "      integer :: rdLocalTime\n");
   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: rdLocal%s\n", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   var_ptr = vars;
   while (var_ptr) {
      fortprintf(fd, "      integer :: rdVarID%s\n", var_ptr->name_in_file);
      var_ptr = var_ptr->next;
   }

   fclose(fd);
   

   /*
    *  Generate read and distribute code
    */
   fd = fopen("io_input_fields.inc", "w");

   var_ptr = vars;
   while (var_ptr) {
      i = 1;
      dimlist_ptr = var_ptr->dimlist;
      if (var_ptr->vtype == INTEGER) sprintf(vtype, "int"); 
      else if (var_ptr->vtype == REAL) sprintf(vtype, "real"); 

      if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
         if (var_ptr->timedim) {
            fortprintf(fd, "      if ((block %% time_levs(1) %% state %% %s %% ioinfo %% input .and. .not. config_do_restart) .or. &\n", var_ptr->super_array);
            fortprintf(fd, "          (block %% time_levs(1) %% state %% %s %% ioinfo %% restart .and. config_do_restart)) then\n", var_ptr->super_array);
         }
         else {
            fortprintf(fd, "      if ((block %% mesh %% %s %% ioinfo %% input .and. .not. config_do_restart) .or. &\n", var_ptr->super_array);
            fortprintf(fd, "          (block %% mesh %% %s %% ioinfo %% restart .and. config_do_restart)) then\n", var_ptr->super_array);
         }
      }
      else {
         if (var_ptr->timedim) {
            fortprintf(fd, "      if ((block %% time_levs(1) %% state %% %s %% ioinfo %% input .and. .not. config_do_restart) .or. &\n", var_ptr->name_in_code);
            fortprintf(fd, "          (block %% time_levs(1) %% state %% %s %% ioinfo %% restart .and. config_do_restart)) then\n", var_ptr->name_in_code);
         }
         else {
            fortprintf(fd, "      if ((block %% mesh %% %s %% ioinfo %% input .and. .not. config_do_restart) .or. &\n", var_ptr->name_in_code);
            fortprintf(fd, "          (block %% mesh %% %s %% ioinfo %% restart .and. config_do_restart)) then\n", var_ptr->name_in_code);
         }
      }
      vert_dim = 0;
      while (dimlist_ptr) {
            if (i < var_ptr->ndims) {
               has_vert_dim = !strcmp( "nVertLevels", dimlist_ptr->dim->name_in_code);
               fortprintf(fd, "      %s%id %% ioinfo %% start(%i) = 1\n", vtype, var_ptr->ndims, i);
               if (has_vert_dim) {
                  vert_dim = i;
                  fortprintf(fd, "#ifdef EXPAND_LEVELS\n");
                  fortprintf(fd, "      if (.not. config_do_restart) then\n");
                  fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = 1\n", vtype, var_ptr->ndims, i);
                  fortprintf(fd, "      else\n");
                  fortprintf(fd, "#endif\n");
               }
               if (dimlist_ptr->dim->constant_value < 0)
                  if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = block %% mesh %% %s\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_code);
                  else fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = block %% mesh %% %s\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_file);
               else
                  fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = %s\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_code);
               if (has_vert_dim) {
                  fortprintf(fd, "#ifdef EXPAND_LEVELS\n");
                  fortprintf(fd, "      end if\n");
                  fortprintf(fd, "#endif\n");
               }
            }
            else {
               if (is_derived_dim(dimlist_ptr->dim->name_in_code)) {
                  split_derived_dim_string(dimlist_ptr->dim->name_in_code, &cp1, &cp2);
                  fortprintf(fd, "      %s%id %% ioinfo %% start(%i) = read%sStart\n", vtype, var_ptr->ndims, i, cp1);
                  fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = read%sCount%s\n", vtype, var_ptr->ndims, i, cp1, cp2);
                  free(cp1);
                  free(cp2);
               }
               else {
                  fortprintf(fd, "      %s%id %% ioinfo %% start(%i) = read%sStart\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_code+1);
                  fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = read%sCount\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_code+1);
               }
            }
         dimlist_ptr = dimlist_ptr->next;
         i++;
      }

      if (var_ptr->ndims > 0) {
         fortprintf(fd, "      allocate(%s%id %% array(", vtype, var_ptr->ndims);
         i = 1;
         dimlist_ptr = var_ptr->dimlist;
   
         if (i < var_ptr->ndims) {
            if (dimlist_ptr->dim->constant_value < 0)
               if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "block %% mesh %% %s", dimlist_ptr->dim->name_in_code);
               else fortprintf(fd, "block %% mesh %% %s", dimlist_ptr->dim->name_in_file);
            else
               fortprintf(fd, "%s", dimlist_ptr->dim->name_in_code);
         }
         else {
            if (is_derived_dim(dimlist_ptr->dim->name_in_code)) {
               split_derived_dim_string(dimlist_ptr->dim->name_in_code, &cp1, &cp2);
               fortprintf(fd, "read%sCount%s", cp1, cp2);
               free(cp1);
               free(cp2);
            }
            else
               fortprintf(fd, "read%sCount", dimlist_ptr->dim->name_in_code+1);
         }
    
         dimlist_ptr = dimlist_ptr->next;
         i++;
         while (dimlist_ptr) {
            if (i < var_ptr->ndims) {
               if (dimlist_ptr->dim->constant_value < 0)
                  if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", block %% mesh %% %s", dimlist_ptr->dim->name_in_code);
                  else fortprintf(fd, ", block %% mesh %% %s", dimlist_ptr->dim->name_in_file);
               else
                  fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_code);
            }
            else {
               if (is_derived_dim(dimlist_ptr->dim->name_in_code)) {
                  split_derived_dim_string(dimlist_ptr->dim->name_in_code, &cp1, &cp2);
                  fortprintf(fd, ", read%sCount%s", cp1, cp2);
                  free(cp1);
                  free(cp2);
               }
               else
                  fortprintf(fd, ", read%sCount", dimlist_ptr->dim->name_in_code+1);
            }
            dimlist_ptr = dimlist_ptr->next;
            i++;
         }
         fortprintf(fd, "))\n\n");

         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            fortprintf(fd, "      allocate(super_%s%id(", vtype, var_ptr->ndims);
            i = 1;
            dimlist_ptr = var_ptr->dimlist;
      
            if (i < var_ptr->ndims) {
               if (dimlist_ptr->dim->constant_value < 0)
                  if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "block %% mesh %% %s", dimlist_ptr->dim->name_in_code);
                  else fortprintf(fd, "block %% mesh %% %s", dimlist_ptr->dim->name_in_file);
               else
                  fortprintf(fd, "%s", dimlist_ptr->dim->name_in_code);
            }
            dimlist_ptr = dimlist_ptr->next;
            i++;
            while (dimlist_ptr) {
               if (dimlist_ptr->dim->constant_value < 0)
                  if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", block %% mesh %% %s", dimlist_ptr->dim->name_in_code);
                  else fortprintf(fd, ", block %% mesh %% %s", dimlist_ptr->dim->name_in_file);
               else
                  fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_code);
               dimlist_ptr = dimlist_ptr->next;
               i++;
            }
            fortprintf(fd, "))\n\n");
         }
      }

      fortprintf(fd, "      %s%id %% ioinfo %% fieldName = \'%s\'\n", vtype, var_ptr->ndims, var_ptr->name_in_file);
      if (var_ptr->timedim)
         fortprintf(fd, "      call io_input_field_time(input_obj, %s%id)\n", vtype, var_ptr->ndims);
      else
         fortprintf(fd, "      call io_input_field(input_obj, %s%id)\n", vtype, var_ptr->ndims);

      if (vert_dim > 0) {
         fortprintf(fd, "#ifdef EXPAND_LEVELS\n");
         fortprintf(fd, "      if (.not. config_do_restart) then\n");
         fortprintf(fd, "         do k=2,EXPAND_LEVELS\n");
         fortprintf(fd, "            %s%id %% array(", vtype, var_ptr->ndims);
         for (i=1; i<=var_ptr->ndims; i++) {
            if (i > 1) fortprintf(fd, ",");
            fortprintf(fd, "%s", i == vert_dim ? "k" : ":");
         }
         fortprintf(fd, ") = %s%id %% array(", vtype, var_ptr->ndims);
         for (i=1; i<=var_ptr->ndims; i++) {
            if (i > 1) fortprintf(fd, ",");
            fortprintf(fd, "%s", i == vert_dim ? "1" : ":");
         }
         fortprintf(fd, ")\n");
         fortprintf(fd, "         end do\n");
         fortprintf(fd, "      end if\n");
         fortprintf(fd, "#endif\n");
      }

      if (var_ptr->ndims > 0) {
         fortprintf(fd, "      call dmpar_alltoall_field(dminfo, &\n");
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            if (var_ptr->timedim) 
               fortprintf(fd, "                                %s%id %% array, super_%s%id, &\n", vtype, var_ptr->ndims, vtype, var_ptr->ndims);
            else
               fortprintf(fd, "                                %s%id %% array, super_%s%id, &\n", vtype, var_ptr->ndims, vtype, var_ptr->ndims);
         }
         else {
            if (var_ptr->timedim) 
               fortprintf(fd, "                                %s%id %% array, block %% time_levs(1) %% state %% %s %% array, &\n", vtype, var_ptr->ndims, var_ptr->name_in_code);
            else
               fortprintf(fd, "                                %s%id %% array, block %% mesh %% %s %% array, &\n", vtype, var_ptr->ndims, var_ptr->name_in_code);
         }
   
         i = 1;
         dimlist_ptr = var_ptr->dimlist;
         
         if (i < var_ptr->ndims)
            if (dimlist_ptr->dim->constant_value < 0)
               if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "                                block %% mesh %% %s", dimlist_ptr->dim->name_in_code);
               else fortprintf(fd, "                                block %% mesh %% %s", dimlist_ptr->dim->name_in_file);
            else
               fortprintf(fd, "                                %s", dimlist_ptr->dim->name_in_code);
         else {
            lastdim = dimlist_ptr;
            if (is_derived_dim(dimlist_ptr->dim->name_in_code)) {
               split_derived_dim_string(dimlist_ptr->dim->name_in_code, &cp1, &cp2);
               fortprintf(fd, "                                read%sCount%s", cp1, cp2);
               free(cp1);
               free(cp2);
            }
            else
               fortprintf(fd, "                                read%sCount", dimlist_ptr->dim->name_in_code+1);
         }
    
         dimlist_ptr = dimlist_ptr->next;
         i++;
         while (dimlist_ptr) {
            if (i < var_ptr->ndims)
               if (dimlist_ptr->dim->constant_value < 0)
                  if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", block %% mesh %% %s", dimlist_ptr->dim->name_in_code);
                  else fortprintf(fd, ", block %% mesh %% %s", dimlist_ptr->dim->name_in_file);
               else
                  fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_code);
            else {
               lastdim = dimlist_ptr;
               if (is_derived_dim(dimlist_ptr->dim->name_in_code)) {
                  split_derived_dim_string(dimlist_ptr->dim->name_in_code, &cp1, &cp2);
                  fortprintf(fd, ", read%sCount%s", cp1, cp2);
                  free(cp1);
                  free(cp2);
               }
               else
                  fortprintf(fd, ", read%sCount", dimlist_ptr->dim->name_in_code+1);
            }
            dimlist_ptr = dimlist_ptr->next;
            i++;
         }
         fortprintf(fd, ", block %% mesh %% %s, &\n", lastdim->dim->name_in_code);
   
         if (is_derived_dim(lastdim->dim->name_in_code)) {
            fortprintf(fd, "                                send%sList, recv%sList)\n", lastdim->dim->name_in_file+1, lastdim->dim->name_in_file+1);
         }
         else
            fortprintf(fd, "                                send%sList, recv%sList)\n", lastdim->dim->name_in_code+1, lastdim->dim->name_in_code+1);


         /* Copy from super_ array to field */
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            i = 1;
            dimlist_ptr = var_ptr->dimlist;
            while (i <= var_ptr->ndims) {
               if (dimlist_ptr->dim->constant_value < 0)
                  if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "      do i%i=1,block %% mesh %% %s\n", i, dimlist_ptr->dim->name_in_code);
                  else fortprintf(fd, "      do i%i=1,block %% mesh %% %s\n", i, dimlist_ptr->dim->name_in_file);
               else
                  fortprintf(fd, "      do i%i=1,%s\n", i, dimlist_ptr->dim->name_in_code);
   
               i++;
               dimlist_ptr = dimlist_ptr->next;
            }
   
            if (var_ptr->timedim) 
               fortprintf(fd, "         block %% time_levs(1) %% state %% %s %% array(index_%s,", var_ptr->super_array, var_ptr->name_in_code);
            else
               fortprintf(fd, "         block %% mesh %% %s %% array(index_%s,", var_ptr->super_array, var_ptr->name_in_code);
            for(i=1; i<=var_ptr->ndims; i++) {
               fortprintf(fd, "i%i",i);
               if (i < var_ptr->ndims) fortprintf(fd, ",");
            }
            fortprintf(fd, ") = super_%s%id(", vtype, var_ptr->ndims);
            for(i=1; i<=var_ptr->ndims; i++) {
               fortprintf(fd, "i%i",i);
               if (i < var_ptr->ndims) fortprintf(fd, ",");
            }
            fortprintf(fd, ")\n");
   
            i = 1;
            while (i <= var_ptr->ndims) {
               fortprintf(fd, "      end do\n");
               i++;
            }
         }

         fortprintf(fd, "      deallocate(%s%id %% array)\n", vtype, var_ptr->ndims);
         if (strncmp(var_ptr->super_array, "-", 1024) != 0)
            fortprintf(fd, "      deallocate(super_%s%id)\n", vtype, var_ptr->ndims);
      }
      else {
         if (var_ptr->timedim) 
            fortprintf(fd, "      block %% time_levs(1) %% state %% %s %% scalar = %s%id %% scalar\n", var_ptr->name_in_code, vtype, var_ptr->ndims);
         else
            fortprintf(fd, "      block %% mesh %% %s %% scalar = %s%id %% scalar\n", var_ptr->name_in_code, vtype, var_ptr->ndims);
      }
     
      fortprintf(fd, "      end if\n\n");

      var_ptr = var_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate NetCDF reads of dimension and variable IDs
    */
   fd = fopen("netcdf_read_ids.inc", "w");

   fortprintf(fd, "      nferr = nf_inq_unlimdim(input_obj %% rd_ncid, input_obj %% rdDimIDTime)\n");
   fortprintf(fd, "      nferr = nf_inq_dimlen(input_obj %% rd_ncid, input_obj %% rdDimIDTime, input_obj %% rdLocalTime)\n");
   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
         fortprintf(fd, "      nferr = nf_inq_dimid(input_obj %% rd_ncid, \'%s\', input_obj %% rdDimID%s)\n", dim_ptr->name_in_file, dim_ptr->name_in_file);
         fortprintf(fd, "      nferr = nf_inq_dimlen(input_obj %% rd_ncid, input_obj %% rdDimID%s, input_obj %% rdLocal%s)\n", dim_ptr->name_in_file, dim_ptr->name_in_file);
      }
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   var_ptr = vars;
   while (var_ptr) {
      fortprintf(fd, "      nferr = nf_inq_varid(input_obj %% rd_ncid, \'%s\', input_obj %% rdVarID%s)\n", var_ptr->name_in_file, var_ptr->name_in_file);
      var_ptr = var_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate code to return dimension given its name
    */
   fd = fopen("get_dimension_by_name.inc", "w");

   dim_ptr = dims;
   while (dim_ptr->constant_value >= 0 || is_derived_dim(dim_ptr->name_in_code)) dim_ptr = dim_ptr->next;
   if (!dim_ptr->namelist_defined) {
      fortprintf(fd, "      if (trim(dimname) == \'%s\') then\n", dim_ptr->name_in_code);
      fortprintf(fd, "         dimsize = input_obj %% rdLocal%s\n", dim_ptr->name_in_file);
   }
   else {
      fortprintf(fd, "      if (trim(dimname) == \'%s\') then\n", dim_ptr->name_in_file);
      fortprintf(fd, "         dimsize = %s\n", dim_ptr->name_in_code);
   }
   dim_ptr = dim_ptr->next;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !is_derived_dim(dim_ptr->name_in_code)) {
         if (!dim_ptr->namelist_defined) {
            fortprintf(fd, "      else if (trim(dimname) == \'%s\') then\n", dim_ptr->name_in_code);
            fortprintf(fd, "         dimsize = input_obj %% rdLocal%s\n", dim_ptr->name_in_file);
         }
         else {
            fortprintf(fd, "      else if (trim(dimname) == \'%s\') then\n", dim_ptr->name_in_file);
            fortprintf(fd, "         dimsize = %s\n", dim_ptr->name_in_code);
         }
      }
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "      end if\n");

   fclose(fd);
   
   
   /*
    *  Generate code to read 0d, 1d, 2d, 3d time-invariant fields
    */
   for(j=0; j<2; j++) {
      for(i=0; i<=3; i++) {
         if (j == 0) {
            sprintf(fname, "input_field%idinteger.inc", i);
            ivtype = INTEGER;
         }
         else {
            sprintf(fname, "input_field%idreal.inc", i);
            ivtype = REAL;
         }
         fd = fopen(fname, "w");
   
         var_ptr = vars;
         while (var_ptr && (var_ptr->ndims != i || var_ptr->vtype != ivtype || var_ptr->timedim)) var_ptr = var_ptr->next;
         if (var_ptr) {
            fortprintf(fd, "      if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
            fortprintf(fd, "         varID = input_obj %% rdVarID%s\n", var_ptr->name_in_file);
            var_ptr = var_ptr->next;
            while (var_ptr) {
               if (var_ptr->ndims == i && var_ptr->vtype == ivtype && !var_ptr->timedim) {
                  fortprintf(fd, "      else if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
                  fortprintf(fd, "         varID = input_obj %% rdVarID%s\n", var_ptr->name_in_file);
               }
               var_ptr = var_ptr->next;
            }
            fortprintf(fd, "      end if\n");
         }
      
         fclose(fd);
      } 
   } 
   
   
   /*
    *  Generate code to read 0d, 1d, 2d, 3d time-varying real fields
    */
   for(i=0; i<=3; i++) { 
      sprintf(fname, "input_field%idreal_time.inc", i);
      fd = fopen(fname, "w");
   
      var_ptr = vars;
      while (var_ptr && (var_ptr->ndims != i || var_ptr->vtype != REAL || !var_ptr->timedim)) var_ptr = var_ptr->next;
      if (var_ptr) {
         fortprintf(fd, "      if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
         fortprintf(fd, "         varID = input_obj %% rdVarID%s\n", var_ptr->name_in_file);
         var_ptr = var_ptr->next;
         while (var_ptr) {
            if (var_ptr->ndims == i && var_ptr->vtype == REAL && var_ptr->timedim) {
               fortprintf(fd, "      else if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
               fortprintf(fd, "         varID = input_obj %% rdVarID%s\n", var_ptr->name_in_file);
            }
            var_ptr = var_ptr->next;
         }
         fortprintf(fd, "      end if\n");
      }
   
      fclose(fd);
   } 
   
}


void gen_writes(struct variable * vars, struct dimension * dims, struct namelist * namelists)
{
   struct variable * var_ptr;
   struct dimension * dim_ptr;
   struct dimension_list * dimlist_ptr, * lastdim;
   struct dtable * dictionary;
   struct namelist * nl;
   FILE * fd;
   char vtype[5];
   char fname[32];
   char * cp1, * cp2;
   int i, j;
   int ivtype;
   
   
   /*
    *  Generate declarations of IDs belonging in io_output_object
    */
   fd = fopen("io_output_obj_decls.inc", "w");

   fortprintf(fd, "      integer :: wrDimIDTime\n");
   dim_ptr = dims;
   while (dim_ptr) {
      fortprintf(fd, "      integer :: wrDimID%s\n", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   var_ptr = vars;
   while (var_ptr) {
      fortprintf(fd, "      integer :: wrVarID%s\n", var_ptr->name_in_file);
      var_ptr = var_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate declarations of temporary dimension variables used for arguments
    */
   fd = fopen("output_dim_actual_decls.inc", "w");

   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: %sGlobal\n", dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: %sGlobal\n", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   fclose(fd);


   /*
    *  Generate initialization of temporary dimension variables used for arguments
    */
   fd = fopen("output_dim_inits.inc", "w");

   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      %sGlobal = block_ptr %% mesh %% %s\n", dim_ptr->name_in_code, dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      %sGlobal = block_ptr %% mesh %% %s\n", dim_ptr->name_in_file, dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   fclose(fd);


   /*
    *  Generate actual dimension argument list
    */
   fd = fopen("output_dim_actual_args.inc", "w");
   dim_ptr = dims;
   if (dim_ptr && dim_ptr->constant_value < 0 && !is_derived_dim(dim_ptr->name_in_code)) {
      if (!dim_ptr->namelist_defined) fortprintf(fd, "                            %sGlobal", dim_ptr->name_in_code);
      else fortprintf(fd, "                            %sGlobal", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, ", %sGlobal", dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, ", %sGlobal", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, " &\n");

   fclose(fd);


   /*
    *  Generate NetCDF calls to define dimensions, variables, and global attributes
    */
   fd = fopen("netcdf_def_dims_vars.inc", "w");

   fortprintf(fd, "      nferr = nf_def_dim(output_obj %% wr_ncid, \'Time\', NF_UNLIMITED, output_obj %% wrDimIDTime)\n");
   dim_ptr = dims;
   while (dim_ptr) {
      fortprintf(fd, "      nferr = nf_def_dim(output_obj %% wr_ncid, \'%s\', %s, output_obj %% wrDimID%s)\n", dim_ptr->name_in_file, dim_ptr->name_in_code, dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   var_ptr = vars;
   while (var_ptr) {
      fortprintf(fd, "      if (.false. &\n");
      if (var_ptr->iostreams & RESTART0) fortprintf(fd, "          .or. output_obj %% stream == RESTART &\n");
      if (var_ptr->iostreams & OUTPUT0)  fortprintf(fd, "          .or. output_obj %% stream == OUTPUT &\n");
      fortprintf(fd, "      ) then\n");
      dimlist_ptr = var_ptr->dimlist;
      i = 1;
      while(dimlist_ptr) {
         fortprintf(fd, "      dimlist(%i) = output_obj %% wrDimID%s\n", i++, dimlist_ptr->dim->name_in_file);
         dimlist_ptr = dimlist_ptr->next;
      }
      if (var_ptr->timedim) fortprintf(fd, "      dimlist(%i) = output_obj %% wrDimIDTime\n", i++);
      if (var_ptr->vtype == INTEGER)
         fortprintf(fd, "      nferr = nf_def_var(output_obj %% wr_ncid, \'%s\', NF_INT, %i, dimlist, output_obj %% wrVarID%s)\n", var_ptr->name_in_file, var_ptr->ndims + var_ptr->timedim, var_ptr->name_in_file);
      else if (var_ptr->vtype == REAL)
         fortprintf(fd, "      nferr = nf_def_var(output_obj %% wr_ncid, \'%s\', NF_DOUBLE, %i, dimlist, output_obj %% wrVarID%s)\n", var_ptr->name_in_file, var_ptr->ndims + var_ptr->timedim, var_ptr->name_in_file);

      fortprintf(fd, "      end if\n\n");

      var_ptr = var_ptr->next;
   }

   nl = namelists;
   while (nl) {
      if (nl->vtype == INTEGER)
         fortprintf(fd, "      nferr = nf_put_att_int(output_obj %% wr_ncid, NF_GLOBAL, \'%s\', NF_INT, 1, %s)\n", nl->name, nl->name);
      else if (nl->vtype == REAL) {
         fortprintf(fd, "      if (RKIND == 8) then\n", nl->name);
         fortprintf(fd, "         nferr = nf_put_att_double(output_obj %% wr_ncid, NF_GLOBAL, \'%s\', NF_DOUBLE, 1, %s)\n", nl->name, nl->name);
         fortprintf(fd, "      else if (RKIND == 4) then\n", nl->name);
         fortprintf(fd, "         nferr = nf_put_att_real(output_obj %% wr_ncid, NF_GLOBAL, \'%s\', NF_FLOAT, 1, %s)\n", nl->name, nl->name);
         fortprintf(fd, "      end if\n");
      }
      else if (nl->vtype == CHARACTER)
         fortprintf(fd, "      nferr = nf_put_att_text(output_obj %% wr_ncid, NF_GLOBAL, \'%s\', len_trim(%s), trim(%s))\n", nl->name, nl->name, nl->name);
      else if (nl->vtype == LOGICAL) {
         fortprintf(fd, "      if (%s) then\n", nl->name);
         fortprintf(fd, "         nferr = nf_put_att_text(output_obj %% wr_ncid, NF_GLOBAL, \'%s\', 1, \'T\')\n", nl->name);
         fortprintf(fd, "      else\n");
         fortprintf(fd, "         nferr = nf_put_att_text(output_obj %% wr_ncid, NF_GLOBAL, \'%s\', 1, \'F\')\n", nl->name);
         fortprintf(fd, "      end if\n");
      }
      nl = nl->next;
   }

   fclose(fd);   
   
   
   /*
    *  Generate collect and write code
    */
   fd = fopen("io_output_fields.inc", "w");

   var_ptr = vars;
   while (var_ptr) {
      i = 1;
      dimlist_ptr = var_ptr->dimlist;
      if (var_ptr->vtype == INTEGER) sprintf(vtype, "int"); 
      else if (var_ptr->vtype == REAL) sprintf(vtype, "real"); 

      if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
         if (var_ptr->timedim) {
            fortprintf(fd, "      if ((domain %% blocklist %% time_levs(1) %% state %% %s %% ioinfo %% output .and. output_obj %% stream == OUTPUT) .or. &\n", var_ptr->super_array);
            fortprintf(fd, "          (domain %% blocklist %% time_levs(1) %% state %% %s %% ioinfo %% restart .and. output_obj %% stream == RESTART)) then\n", var_ptr->super_array);
         }
         else {
            fortprintf(fd, "      if ((domain %% blocklist %% mesh %% %s %% ioinfo %% output .and. output_obj %% stream == OUTPUT) .or. &\n", var_ptr->super_array);
            fortprintf(fd, "          (domain %% blocklist %% mesh %% %s %% ioinfo %% restart .and. output_obj %% stream == RESTART)) then\n", var_ptr->super_array);
         }
      }
      else {
         if (var_ptr->timedim) {
            fortprintf(fd, "      if ((domain %% blocklist %% time_levs(1) %% state %% %s %% ioinfo %% output .and. output_obj %% stream == OUTPUT) .or. &\n", var_ptr->name_in_code);
            fortprintf(fd, "          (domain %% blocklist %% time_levs(1) %% state %% %s %% ioinfo %% restart .and. output_obj %% stream == RESTART)) then\n", var_ptr->name_in_code);
         }
         else {
            fortprintf(fd, "      if ((domain %% blocklist %% mesh %% %s %% ioinfo %% output .and. output_obj %% stream == OUTPUT) .or. &\n", var_ptr->name_in_code);
            fortprintf(fd, "          (domain %% blocklist %% mesh %% %s %% ioinfo %% restart .and. output_obj %% stream == RESTART)) then\n", var_ptr->name_in_code);
         }
      }

      if (var_ptr->ndims > 0) {
         while (dimlist_ptr) {
               if (i < var_ptr->ndims) {
                  fortprintf(fd, "      %s%id %% ioinfo %% start(%i) = 1\n", vtype, var_ptr->ndims, i);
                  if (dimlist_ptr->dim->constant_value < 0)
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = domain %% blocklist %% mesh %% %s\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_code);
                     else fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = domain %% blocklist %% mesh %% %s\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_file);
                  else
                     fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = %s\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_code);
               }
               else {
                  fortprintf(fd, "      %s%id %% ioinfo %% start(%i) = 1\n", vtype, var_ptr->ndims, i);
                  if (is_derived_dim(dimlist_ptr->dim->name_in_code)) {
                     split_derived_dim_string(dimlist_ptr->dim->name_in_code, &cp1, &cp2);
                     fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = n%sGlobal%s\n", vtype, var_ptr->ndims, i, cp1, cp2);
                     free(cp1);
                     free(cp2);
                  }
                  else
                     fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = %sGlobal\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_code);
               }
            dimlist_ptr = dimlist_ptr->next;
            i++;
         }
   
         fortprintf(fd, "      allocate(%s%id %% array(", vtype, var_ptr->ndims);
         i = 1;
         dimlist_ptr = var_ptr->dimlist;
   
         if (i < var_ptr->ndims)
            if (dimlist_ptr->dim->constant_value < 0)
               if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_code);
               else fortprintf(fd, "domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_file);
            else
               fortprintf(fd, "%s", dimlist_ptr->dim->name_in_code);
         else {
            if (is_derived_dim(dimlist_ptr->dim->name_in_code)) {
               split_derived_dim_string(dimlist_ptr->dim->name_in_code, &cp1, &cp2);
               fortprintf(fd, "n%sGlobal%s", cp1, cp2);
               free(cp1);
               free(cp2);
            }
            else
               fortprintf(fd, "%sGlobal", dimlist_ptr->dim->name_in_code);
            lastdim = dimlist_ptr;
         }
         dimlist_ptr = dimlist_ptr->next;
         i++;
         while (dimlist_ptr) {
            if (i < var_ptr->ndims)
               if (dimlist_ptr->dim->constant_value < 0)
                  if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_code);
                  else fortprintf(fd, ", domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_file);
               else
                  fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_code);
            else {
               if (is_derived_dim(dimlist_ptr->dim->name_in_code)) {
                  split_derived_dim_string(dimlist_ptr->dim->name_in_code, &cp1, &cp2);
                  fortprintf(fd, ", n%sGlobal%s", cp1, cp2);
                  free(cp1);
                  free(cp2);
               }
               else
                  fortprintf(fd, ", %sGlobal", dimlist_ptr->dim->name_in_code);
               lastdim = dimlist_ptr;
            }
            dimlist_ptr = dimlist_ptr->next;
            i++;
         }
         fortprintf(fd, "))\n\n");

         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            if (var_ptr->ndims > 0) {
               fortprintf(fd, "      allocate(super_%s%id(", vtype, var_ptr->ndims);
               i = 1;
               dimlist_ptr = var_ptr->dimlist;
               while (dimlist_ptr) {
                  if (dimlist_ptr->dim->constant_value < 0)
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_code);
                     else fortprintf(fd, "domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_file);
                  else
                     fortprintf(fd, "%s", dimlist_ptr->dim->name_in_code);
   
                  if (i < var_ptr->ndims) fortprintf(fd, ", ");
      
                  dimlist_ptr = dimlist_ptr->next;
                  i++;
               }
               fortprintf(fd, "))\n\n");
            }

            /* Copy from field to super_ array */
            i = 1;
            dimlist_ptr = var_ptr->dimlist;
            while (i <= var_ptr->ndims) {
               if (dimlist_ptr->dim->constant_value < 0)
                  if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "      do i%i=1,domain %% blocklist %% mesh %% %s\n", i, dimlist_ptr->dim->name_in_code);
                  else fortprintf(fd, "      do i%i=1,domain %% blocklist %% mesh %% %s\n", i, dimlist_ptr->dim->name_in_file);
               else
                  fortprintf(fd, "      do i%i=1,%s\n", i, dimlist_ptr->dim->name_in_code);

               i++;
               dimlist_ptr = dimlist_ptr->next;
            }

            fortprintf(fd, "         super_%s%id(", vtype, var_ptr->ndims);
            for(i=1; i<=var_ptr->ndims; i++) {
               fortprintf(fd, "i%i",i);
               if (i < var_ptr->ndims) fortprintf(fd, ",");
            }
            if (var_ptr->timedim) 
               fortprintf(fd, ") = domain %% blocklist %% time_levs(1) %% state %% %s %% array(", var_ptr->super_array);
            else
               fortprintf(fd, ") = domain %% blocklist %% mesh %% %s %% array(", var_ptr->super_array);
            fortprintf(fd, "index_%s", var_ptr->name_in_code);
            for(i=1; i<=var_ptr->ndims; i++) {
               fortprintf(fd, ",i%i",i);
            }
            fortprintf(fd, ")\n");

            i = 1;
            while (i <= var_ptr->ndims) {
               fortprintf(fd, "      end do\n");
               i++;
            }
         }

         fortprintf(fd, "      %s%id %% ioinfo %% fieldName = \'%s\'\n", vtype, var_ptr->ndims, var_ptr->name_in_file);
         fortprintf(fd, "      call dmpar_alltoall_field(domain %% dminfo, &\n");
         if (strncmp(var_ptr->super_array, "-", 1024) != 0)
            fortprintf(fd, "                                super_%s%id, %s%id %% array, &\n", vtype, var_ptr->ndims, vtype, var_ptr->ndims);
         else {
            if (var_ptr->timedim) 
               fortprintf(fd, "                                domain %% blocklist %% time_levs(1) %% state %% %s %% array, %s%id %% array, &\n", var_ptr->name_in_code, vtype, var_ptr->ndims);
            else
               fortprintf(fd, "                                domain %% blocklist %% mesh %% %s %% array, %s%id %% array, &\n", var_ptr->name_in_code, vtype, var_ptr->ndims);
         }
   
         i = 1;
         dimlist_ptr = var_ptr->dimlist;
         
         if (dimlist_ptr->dim->constant_value < 0)
            if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "                                domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_code);
            else fortprintf(fd, "                                domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_file);
         else
            fortprintf(fd, "                                %s", dimlist_ptr->dim->name_in_code);
    
         dimlist_ptr = dimlist_ptr->next;
         i++;
         while (dimlist_ptr) {
            if (dimlist_ptr->dim->constant_value < 0)
               if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_code);
               else fortprintf(fd, ", domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_file);
            else
               fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_code);
   
            dimlist_ptr = dimlist_ptr->next;
            i++;
         }     
   
         if (is_derived_dim(lastdim->dim->name_in_code)) {
            split_derived_dim_string(lastdim->dim->name_in_code, &cp1, &cp2);
            fortprintf(fd, ", n%sGlobal%s, &\n", cp1, cp2);
            fortprintf(fd, "                                output_obj %% send%sList, output_obj %% recv%sList)\n", lastdim->dim->name_in_file+1, lastdim->dim->name_in_file+1);
            free(cp1);
            free(cp2);
         }
         else {
            fortprintf(fd, ", %sGlobal, &\n", lastdim->dim->name_in_code);
            fortprintf(fd, "                                output_obj %% send%sList, output_obj %% recv%sList)\n", lastdim->dim->name_in_code+1, lastdim->dim->name_in_code+1);
         }
      }
      else {
         fortprintf(fd, "      %s%id %% ioinfo %% fieldName = \'%s\'\n", vtype, var_ptr->ndims, var_ptr->name_in_file);
         if (var_ptr->timedim) 
            fortprintf(fd, "      %s%id %% scalar = domain %% blocklist %% time_levs(1) %% state %% %s %% scalar\n", vtype, var_ptr->ndims, var_ptr->name_in_code);
         else
            fortprintf(fd, "      %s%id %% scalar = domain %% blocklist %% mesh %% %s %% scalar\n", vtype, var_ptr->ndims, var_ptr->name_in_code);
      }

      if (var_ptr->timedim)
         fortprintf(fd, "      if (domain %% dminfo %% my_proc_id == IO_NODE) call io_output_field_time(output_obj, %s%id)\n", vtype, var_ptr->ndims);
      else
         fortprintf(fd, "      if (domain %% dminfo %% my_proc_id == IO_NODE) call io_output_field(output_obj, %s%id)\n", vtype, var_ptr->ndims);
      if (var_ptr->ndims > 0) {
         fortprintf(fd, "      deallocate(%s%id %% array)\n", vtype, var_ptr->ndims);
         if (strncmp(var_ptr->super_array, "-", 1024) != 0)
            fortprintf(fd, "      deallocate(super_%s%id)\n", vtype, var_ptr->ndims);
      }
      fortprintf(fd, "      end if\n\n");

      var_ptr = var_ptr->next;
   }

   fclose(fd);
   
   
   /*
    *  Generate code to write 0d, 1d, 2d, 3d time-invariant fields
    */
   for(j=0; j<2; j++) {
      for(i=0; i<=3; i++) {
         if (j == 0) {
            sprintf(fname, "output_field%idinteger.inc", i);
            ivtype = INTEGER;
         }
         else {
            sprintf(fname, "output_field%idreal.inc", i);
            ivtype = REAL;
         }
         fd = fopen(fname, "w");
   
         var_ptr = vars;
         while (var_ptr && (var_ptr->ndims != i || var_ptr->vtype != ivtype || var_ptr->timedim)) var_ptr = var_ptr->next;
         if (var_ptr) {
            fortprintf(fd, "      if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
            fortprintf(fd, "         varID = output_obj %% wrVarID%s\n", var_ptr->name_in_file);
            var_ptr = var_ptr->next;
            while (var_ptr) {
               if (var_ptr->ndims == i && var_ptr->vtype == ivtype && !var_ptr->timedim) {
                  fortprintf(fd, "      else if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
                  fortprintf(fd, "         varID = output_obj %% wrVarID%s\n", var_ptr->name_in_file);
               }
               var_ptr = var_ptr->next;
            }
            fortprintf(fd, "      end if\n");
         }
      
         fclose(fd);
      } 
   } 

   
   /*
    *  Generate code to write 0d, 1d, 2d, 3d real time-varying fields
    */
   for(i=0; i<=3; i++) {
      sprintf(fname, "output_field%idreal_time.inc", i);
      fd = fopen(fname, "w");

      var_ptr = vars;
      while (var_ptr && (var_ptr->ndims != i || var_ptr->vtype != REAL || !var_ptr->timedim)) var_ptr = var_ptr->next;
      if (var_ptr) {
         fortprintf(fd, "      if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
         fortprintf(fd, "         varID = output_obj %% wrVarID%s\n", var_ptr->name_in_file);
         var_ptr = var_ptr->next;
         while (var_ptr) {
            if (var_ptr->ndims == i && var_ptr->vtype == REAL && var_ptr->timedim) {
               fortprintf(fd, "      else if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
               fortprintf(fd, "         varID = output_obj %% wrVarID%s\n", var_ptr->name_in_file);
            }
            var_ptr = var_ptr->next;
         }
         fortprintf(fd, "      end if\n");
      }
   
      fclose(fd);
   }
   
}
