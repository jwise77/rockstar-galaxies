#ifndef _IO_HDF5_H_
#define _IO_HDF5_H_
#ifdef ENABLE_HDF5
#include <hdf5.h> /* HDF5 required */
#include <inttypes.h>

hid_t check_H5Fopen(char *filename, unsigned flags);

hid_t check_H5Gopen(hid_t HDF_FileID, char *gid, char *filename);

hid_t check_H5Dopen(hid_t HDF_GroupID, char *dataid, char *gid, char *filename);
hid_t check_H5Dget_space(hid_t HDF_DatasetID);
void check_H5Dread(hid_t HDF_DatasetID, hid_t type, void *buffer, char *dataid, char *gid, char *filename);

hid_t check_H5Aopen_name(hid_t HDF_GroupID, char *dataid, char *gid, char *filename);
hid_t check_H5Aget_space(hid_t HDF_AttrID);
void check_H5Aread(hid_t HDF_AttrID, hid_t type, void *buffer, char *dataid, char *gid, char *filename);

void check_H5Sselect_all(hid_t HDF_DataspaceID);
int64_t check_H5Sget_simple_extent_ndims(hid_t HDF_DataspaceID);
void check_H5Sget_simple_extent_dims(hid_t HDF_DataspaceID, hsize_t *dimsize);

#endif /* ENABLE_HDF5 */
#endif /* _IO_HDF5_H_ */

