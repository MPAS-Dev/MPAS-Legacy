#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "registry_types.h"
#include "gen_inc.h"

int parse_reg(FILE *, struct namelist **, struct dimension **, struct variable **, struct group_list **);
int getword(FILE *, char *);
int is_integer_constant(char *);
void sort_vars(struct variable *);
void sort_group_vars(struct group_list *);

int main(int argc, char ** argv)
{
   FILE * regfile;
   struct namelist * nls;
   struct dimension * dims;
   struct variable * vars;
   struct group_list * groups;

   if (argc != 2) {
      fprintf(stderr,"\nUsage: %s filename\n\n", argv[0]);
      return 1;
   }

   if (regfile = fopen(argv[1], "r")) {
      nls = NULL;
      dims = NULL;
      vars = NULL;
      if (parse_reg(regfile, &nls, &dims, &vars, &groups)) {
         return 1;
      }
   }   
   else {
      fprintf(stderr,"\nError: Could not open file %s for reading.\n\n", argv[1]);
      return 1;
   }   

   sort_vars(vars);
   sort_group_vars(groups);

   gen_namelists(nls);
   gen_field_defs(groups, vars, dims);
   gen_reads(groups, vars, dims);
   gen_writes(groups, vars, dims, nls);

   return 0;
}


int parse_reg(FILE * regfile, struct namelist ** nls, struct dimension ** dims, struct variable ** vars, struct group_list ** groups)
{
   char word[1024];
   struct namelist * nls_ptr;
   struct namelist * nls_chk_ptr;
   struct dimension * dim_ptr;
   struct variable * var_ptr;
   struct dimension_list * dimlist_ptr;
   struct dimension * dimlist_cursor;
   struct group_list * grouplist_ptr;
   struct variable_list * vlist_cursor;

   NEW_NAMELIST(nls_ptr)
   NEW_DIMENSION(dim_ptr)
   NEW_VARIABLE(var_ptr)
   NEW_GROUP_LIST(grouplist_ptr);
   *nls = nls_ptr;
   *dims = dim_ptr;
   *vars = var_ptr;
   *groups = grouplist_ptr;

   while(getword(regfile, word) != EOF) {
      if (strncmp(word, "namelist", 1024) == 0) {
         NEW_NAMELIST(nls_ptr->next)
         nls_ptr = nls_ptr->next;

         getword(regfile, word); 
         if (strncmp(word, "real", 1024) == 0) 
            nls_ptr->vtype = REAL;
         else if (strncmp(word, "integer", 1024) == 0) 
            nls_ptr->vtype = INTEGER;
         else if (strncmp(word, "logical", 1024) == 0) 
            nls_ptr->vtype = LOGICAL;
         else if (strncmp(word, "character", 1024) == 0) 
            nls_ptr->vtype = CHARACTER;

         getword(regfile, nls_ptr->record); 
         getword(regfile, nls_ptr->name); 

         getword(regfile, word); 
         if (nls_ptr->vtype == REAL) 
            nls_ptr->defval.rval = (float)atof(word);
         else if (nls_ptr->vtype == INTEGER) 
            nls_ptr->defval.ival = atoi(word);
         else if (nls_ptr->vtype == LOGICAL) {
            if (strncmp(word, "true", 1024) == 0) 
               nls_ptr->defval.lval = 1;
            else if (strncmp(word, "false", 1024) == 0) 
               nls_ptr->defval.lval = 0;
         }
         else if (nls_ptr->vtype == CHARACTER) 
            strncpy(nls_ptr->defval.cval, word, 32);
      }
      else if (strncmp(word, "dim", 1024) == 0) {
         NEW_DIMENSION(dim_ptr->next)
         dim_ptr = dim_ptr->next;
         dim_ptr->namelist_defined = 0;
         getword(regfile, dim_ptr->name_in_file); 
         getword(regfile, dim_ptr->name_in_code); 
         dim_ptr->constant_value = is_integer_constant(dim_ptr->name_in_code);
         if (strncmp(dim_ptr->name_in_code, "namelist:", 9) == 0) {
            dim_ptr->namelist_defined = 1;
            sprintf(dim_ptr->name_in_code, "%s", (dim_ptr->name_in_code)+9);
            
            /* Check that the referenced namelist variable is defined as an integer variable */
            nls_chk_ptr = (*nls)->next;
            while (nls_chk_ptr) {
               if (strncmp(nls_chk_ptr->name, dim_ptr->name_in_code, 1024) == 0) {
                  if (nls_chk_ptr->vtype != INTEGER) {
                     printf("\nRegistry error: Namelist variable %s must be an integer for namelist-derived dimension %s\n\n", nls_chk_ptr->name, dim_ptr->name_in_file);
                     return 1;
                  }
                  break;
               } 
               nls_chk_ptr = nls_chk_ptr->next;
            }
            if (!nls_chk_ptr) {
               printf("\nRegistry error: Namelist variable %s not defined for namelist-derived dimension %s\n\n", dim_ptr->name_in_code, dim_ptr->name_in_file);
               return 1;
            }
         }
      }
      else if (strncmp(word, "var", 1024) == 0) {
         NEW_VARIABLE(var_ptr->next)
         var_ptr = var_ptr->next;
         var_ptr->ndims = 0;
         var_ptr->timedim = 0;
         var_ptr->iostreams = 0;

         /* 
          * persistence 
          */
         getword(regfile, word); 
         if (strncmp(word, "persistent", 1024) == 0) 
            var_ptr->persistence = PERSISTENT;
         else if (strncmp(word, "scratch", 1024) == 0) 
            var_ptr->persistence = SCRATCH;

         getword(regfile, word); 
         if (strncmp(word, "real", 1024) == 0) 
            var_ptr->vtype = REAL;
         else if (strncmp(word, "integer", 1024) == 0) 
            var_ptr->vtype = INTEGER;
         else if (strncmp(word, "logical", 1024) == 0) 
            var_ptr->vtype = LOGICAL;
         else if (strncmp(word, "text", 1024) == 0) 
            var_ptr->vtype = CHARACTER;

         getword(regfile, var_ptr->name_in_file); 

         NEW_DIMENSION_LIST(dimlist_ptr)
         var_ptr->dimlist = dimlist_ptr;

         getword(regfile, word); /* Should have just read a right paren */
         getword(regfile, word); 
         while (strncmp(word, ")", 1024) != 0) {
            
            if (strncmp(word, "Time", 1024) == 0) {
               var_ptr->timedim = 1;
            }
            else {
               NEW_DIMENSION_LIST(dimlist_ptr->next)
               dimlist_ptr->next->prev = dimlist_ptr;
               dimlist_ptr = dimlist_ptr->next;

               dimlist_cursor = (*dims)->next;
               while (dimlist_cursor && (strncmp(word, dimlist_cursor->name_in_file, 1024) != 0)) dimlist_cursor = dimlist_cursor->next;
               if (dimlist_cursor) {
                  dimlist_ptr->dim = dimlist_cursor;
               }
               else {
                  fprintf(stderr, "Error: Unknown dimension %s for variable %s\n", word, var_ptr->name_in_file);
                  return 1;
               }
            }
            getword(regfile, word); 
         }

         /* 
          * time_dim 
          */
         getword(regfile, word);
         var_ptr->ntime_levs = atoi(word);

         /* 
          * I/O info 
          */
         getword(regfile, word);
         if (strchr(word, (int)'i')) var_ptr->iostreams |= INPUT0;
         if (strchr(word, (int)'r')) var_ptr->iostreams |= RESTART0;
         if (strchr(word, (int)'o')) var_ptr->iostreams |= OUTPUT0;

         getword(regfile, var_ptr->name_in_code); 

         /* 
          * struct 
          */
         getword(regfile, var_ptr->struct_group); 
         grouplist_ptr = *groups;
         grouplist_ptr = grouplist_ptr->next;
         while (grouplist_ptr && strncmp(var_ptr->struct_group, grouplist_ptr->name, 1024)) {
            grouplist_ptr = grouplist_ptr->next;
         }
         if (!grouplist_ptr) {
            grouplist_ptr = *groups;
            while(grouplist_ptr->next) grouplist_ptr = grouplist_ptr->next;
            NEW_GROUP_LIST(grouplist_ptr->next);
            grouplist_ptr = grouplist_ptr->next;
            memcpy(grouplist_ptr->name, var_ptr->struct_group, (size_t)1024);
            NEW_VARIABLE_LIST(grouplist_ptr->vlist);
            grouplist_ptr->vlist->var = var_ptr;
         }
         else {
            vlist_cursor = grouplist_ptr->vlist;
            while (vlist_cursor->next) vlist_cursor = vlist_cursor->next;
            NEW_VARIABLE_LIST(vlist_cursor->next);
            vlist_cursor->next->prev = vlist_cursor;
            vlist_cursor = vlist_cursor->next;
            vlist_cursor->var = var_ptr;
         }


         getword(regfile, var_ptr->super_array);
         getword(regfile, var_ptr->array_class);

         dimlist_ptr = var_ptr->dimlist;
         if (var_ptr->dimlist) var_ptr->dimlist = var_ptr->dimlist->next;
         if (dimlist_ptr) free(dimlist_ptr);

         dimlist_ptr = var_ptr->dimlist;
         while (dimlist_ptr) {
            var_ptr->ndims++; 
            dimlist_ptr = dimlist_ptr->next;
         }
      }
   } 

   nls_ptr = *nls;
   if ((*nls)->next) *nls = (*nls)->next;
   if (nls_ptr) free(nls_ptr);

   dim_ptr = *dims;
   if ((*dims)->next) *dims = (*dims)->next;
   if (dim_ptr) free(dim_ptr);

   var_ptr = *vars;
   if ((*vars)->next) *vars = (*vars)->next;
   if (var_ptr) free(var_ptr);

   grouplist_ptr = *groups;
   if ((*groups)->next) *groups = (*groups)->next;
   if (grouplist_ptr) free(grouplist_ptr);

   return 0;
}

int getword(FILE * regfile, char * word)
{
   int i;
   int c;

   i = 0;
   
   do { c = getc(regfile); } while (((char)c == ' ' || (char)c == '\n' || (char)c == '\t') && c != EOF);

   while ((char)c == '#') {
      do { c = getc(regfile); } while ((char)c != '\n' && c != EOF);
      do { c = getc(regfile); } while (((char)c == ' ' || (char)c == '\n' || (char)c == '\t') && c != EOF);
   };
   while((char)c != ' ' && (char)c != '\n' && (char)c != '\t' && c != EOF && (char)c != '#') {
      word[i++] = (char)c; 
      c = (char)getc(regfile);
   } 
   word[i] = '\0';

   if ((char)c == '#') do { c = getc(regfile); } while ((char)c != '\n' && c != EOF);

   return c;
}

int is_integer_constant(char * c) {
   int i;

   i = 0;
   while (c[i] != '\0') {
      if (c[i] < '0' || c[i] > '9') return -1;
      i++;
   }

   return atoi(c);
}

void sort_vars(struct variable * vars)
{
   struct variable * var_ptr;
   struct variable * var_ptr2;
   struct variable * var_ptr2_prev;
   char super_array[1024];
   char array_class[1024];

   var_ptr = vars;

/* Attempt at sorting first on super-array, then on class in the same loop
   while (var_ptr) {
      memcpy(super_array, var_ptr->super_array, 1024);
      memcpy(array_class, var_ptr->array_class, 1024);
      var_ptr2_prev = var_ptr;
      var_ptr2 = var_ptr->next;
      if (var_ptr2 && 
          (strncmp(super_array, var_ptr2->super_array, 1024) != 0 || strncmp(array_class, var_ptr2->array_class, 1024) != 0)) {
         while (var_ptr2) {
            if (strncmp(super_array, var_ptr2->super_array, 1024) == 0 && strncmp(array_class, var_ptr2->array_class, 1024) == 0) {
               var_ptr2_prev->next = var_ptr2->next;
               var_ptr2->next = var_ptr->next;
               var_ptr->next = var_ptr2;
               var_ptr2 = var_ptr2_prev->next;
            }
            else {
               var_ptr2_prev = var_ptr2_prev->next;
               var_ptr2 = var_ptr2->next;
            }
         }
      } 
      var_ptr = var_ptr->next;
   }
*/

   while (var_ptr) {
      memcpy(super_array, var_ptr->super_array, 1024);
      var_ptr2_prev = var_ptr;
      var_ptr2 = var_ptr->next;
      if (var_ptr2 && strncmp(super_array, var_ptr2->super_array, 1024) != 0) {
         while (var_ptr2) {
            if (strncmp(super_array, var_ptr2->super_array, 1024) == 0) {
               var_ptr2_prev->next = var_ptr2->next;
               var_ptr2->next = var_ptr->next;
               var_ptr->next = var_ptr2;
               var_ptr2 = var_ptr2_prev->next;
            }
            else {
               var_ptr2_prev = var_ptr2_prev->next;
               var_ptr2 = var_ptr2->next;
            }
         }
      } 
      var_ptr = var_ptr->next;
   }

   var_ptr = vars;

   while (var_ptr) {
      memcpy(array_class, var_ptr->array_class, 1024);
      var_ptr2_prev = var_ptr;
      var_ptr2 = var_ptr->next;
      if (var_ptr2 && strncmp(array_class, var_ptr2->array_class, 1024) != 0) {
         while (var_ptr2) {
            if (strncmp(array_class, var_ptr2->array_class, 1024) == 0) {
               var_ptr2_prev->next = var_ptr2->next;
               var_ptr2->next = var_ptr->next;
               var_ptr->next = var_ptr2;
               var_ptr2 = var_ptr2_prev->next;
            }
            else {
               var_ptr2_prev = var_ptr2_prev->next;
               var_ptr2 = var_ptr2->next;
            }
         }
      } 
      var_ptr = var_ptr->next;
   }
}


void sort_group_vars(struct group_list * groups)
{
   struct variable_list * var_list;
   struct variable_list * var_ptr;
   struct variable_list * var_ptr2;
   struct variable_list * var_ptr2_prev;
   struct group_list * group_ptr;
   char super_array[1024];
   char array_class[1024];

   group_ptr = groups;

   while (group_ptr) {

      var_ptr = group_ptr->vlist;
   
      while (var_ptr) {
         memcpy(super_array, var_ptr->var->super_array, 1024);
         var_ptr2_prev = var_ptr;
         var_ptr2 = var_ptr->next;
         if (var_ptr2 && strncmp(super_array, var_ptr2->var->super_array, 1024) != 0) {
            while (var_ptr2) {
               if (strncmp(super_array, var_ptr2->var->super_array, 1024) == 0) {
                  var_ptr2_prev->next = var_ptr2->next;
                  var_ptr2->next = var_ptr->next;
                  var_ptr->next = var_ptr2;
                  var_ptr2 = var_ptr2_prev->next;
               }
               else {
                  var_ptr2_prev = var_ptr2_prev->next;
                  var_ptr2 = var_ptr2->next;
               }
            }
         } 
         var_ptr = var_ptr->next;
      }
   
      var_ptr = group_ptr->vlist;
   
      while (var_ptr) {
         memcpy(array_class, var_ptr->var->array_class, 1024);
         var_ptr2_prev = var_ptr;
         var_ptr2 = var_ptr->next;
         if (var_ptr2 && strncmp(array_class, var_ptr2->var->array_class, 1024) != 0) {
            while (var_ptr2) {
               if (strncmp(array_class, var_ptr2->var->array_class, 1024) == 0) {
                  var_ptr2_prev->next = var_ptr2->next;
                  var_ptr2->next = var_ptr->next;
                  var_ptr->next = var_ptr2;
                  var_ptr2 = var_ptr2_prev->next;
               }
               else {
                  var_ptr2_prev = var_ptr2_prev->next;
                  var_ptr2 = var_ptr2->next;
               }
            }
         } 
         var_ptr = var_ptr->next;
      }

      group_ptr = group_ptr->next;
   }
}
