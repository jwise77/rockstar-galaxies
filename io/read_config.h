#ifndef READ_CONFIG_H
#define READ_CONFIG_H

struct configfile {
  int num_entries;
  char **keys;
  char **values;
  int *touched;
};

char *config_to_string(struct configfile *c, char *key, char *def_val);
double config_to_real(struct configfile *c, char *key, double def_val);
void config_to_real3(struct configfile *c, char *key, double *res, char *def_val);
void syntax_check(struct configfile *c, char *prefix);
void load_config(struct configfile *c, char *filename);
void write_config(struct configfile c, char *filename);
void free_config(struct configfile c);

#endif /* READ_CONFIG_H */
