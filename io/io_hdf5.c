#ifdef ENABLE_HDF5
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h> /* HDF5 required */
#include <inttypes.h>

hid_t check_H5Fopen(char *filename, unsigned flags) {
  hid_t HDF_FileID = H5Fopen(filename, flags, H5P_DEFAULT);
  if (HDF_FileID < 0) {
    fprintf(stderr, "[Error] Failed to open HDF5 file %s!\n", filename);
    exit(1);
  }
  return HDF_FileID;
}

hid_t check_H5Gopen(hid_t HDF_FileID, char *gid, char *filename) {
  hid_t HDF_GroupID = H5Gopen(HDF_FileID, gid);
  if (HDF_GroupID < 0) {
    fprintf(stderr, "[Error] Failed to open group %s in HDF5 file %s!\n", gid, filename);
    exit(1);
  }
  return HDF_GroupID;
}

hid_t check_H5Dopen(hid_t HDF_GroupID, char *dataid, char *gid, char *filename){
  hid_t HDF_DatasetID = H5Dopen(HDF_GroupID, dataid);
  if (HDF_DatasetID < 0) {
    fprintf(stderr, "[Error] Failed to open dataset %s/%s in HDF5 file %s!\n", gid, dataid, filename);
    exit(1);
  }
  return HDF_DatasetID;
}

hid_t check_H5Dget_space(hid_t HDF_DatasetID) {
  hid_t HDF_DataspaceID = H5Dget_space(HDF_DatasetID);
  if (HDF_DataspaceID < 0) {
    fprintf(stderr, "[Error] Failed to get HDF5 dataspace!\n");
    exit(1);
  }
  return HDF_DataspaceID;
}

void check_H5Dread(hid_t HDF_DatasetID, hid_t type, void *buffer, char *dataid, char *gid, char *filename) {
  if (H5Dread(HDF_DatasetID, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer) < 0){
    fprintf(stderr, "[Error] failed to read dataspace %s/%s in HDF5 file %s\n", gid, dataid, filename);
    exit(1);
  }
}

hid_t check_H5Aopen_name(hid_t HDF_GroupID, char *dataid, char *gid, char *filename){
  hid_t HDF_AttrID = H5Aopen_name(HDF_GroupID, dataid);
  if (HDF_AttrID < 0) {
    fprintf(stderr, "[Error] Failed to open attribute %s/%s in HDF5 file %s!\n", gid, dataid, filename);
    exit(1);
  }
  return HDF_AttrID;
}

hid_t check_H5Aget_space(hid_t HDF_AttrID) {
  hid_t HDF_DataspaceID = H5Aget_space(HDF_AttrID);
  if (HDF_AttrID < 0) {
    fprintf(stderr, "[Error] Failed to get HDF5 dataspace!\n");
    exit(1);
  }
  return HDF_DataspaceID;
}

void check_H5Aread(hid_t HDF_AttrID, hid_t type, void *buffer, char *dataid, char *gid, char *filename) {
  if (H5Aread(HDF_AttrID, type, buffer) < 0) {
    fprintf(stderr, "[Error] failed to read attribute %s/%s in HDF5 file %s\n", gid, dataid, filename);
    exit(1);
  }
}

void check_H5Sselect_all(hid_t HDF_DataspaceID) {
  if (H5Sselect_all(HDF_DataspaceID)<0 || H5Sselect_valid(HDF_DataspaceID)<=0) {
    fprintf(stderr, "[Error] Failed to select all elements in HDF5 dataspace!\n");
    exit(1);
  }
}

int64_t check_H5Sget_simple_extent_ndims(hid_t HDF_DataspaceID) {
  int64_t ndims = H5Sget_simple_extent_ndims(HDF_DataspaceID);
  if (ndims < 0) {
    fprintf(stderr, "[Error] Failed to get number of dimensions for HDF5 dataspace!\n");
    exit(1);
  }
  return ndims;
}

void check_H5Sget_simple_extent_dims(hid_t HDF_DataspaceID, hsize_t *dimsize) {
  if (H5Sget_simple_extent_dims(HDF_DataspaceID, dimsize, NULL)<0) {
    fprintf(stderr, "[Error] Failed to get dimensions for HDF5 dataspace!\n");
    exit(1);
  }
}

#endif /* ENABLE_HDF5 */
