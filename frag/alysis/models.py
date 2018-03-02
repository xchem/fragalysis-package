import uuid


def split_res_water():
    pass
class StructHolder(object):
    """
    A class to holad
    """
    def __init__(self, struct_id, resid_pdb=None, ligand=None, water_pdb=None, apo_pdb=None):
        self.struct_id = struct_id
        self.water_pdb = water_pdb
        self.resid_pdb = resid_pdb
        self.ligand = ligand
        if resid_pdb or water_pdb == None:
            new_water_pdb, new_resid_pdb = split_res_water()
            if resid_pdb == None:
                self.resid_pdb == new_resid_pdb
            if water_pdb == None:
                self.water_pdb = new_water_pdb

class ClusterStuff(object):
    """
    A class to cluster stuff.
    Probably also need a write output function too. Once I've got my head around that...
    """
    def __init__(self, parser,lamb,cluster):
        self.parser = parser
        self.lamb = lamb
        self.cluster = cluster
        self.owner_list = []
        self.out_clusters = {}

    def run(self,input_data):
        # self.parser returns a list of Owner objects -
        self.owner_list = self.parser(input_data)
        # Now convert o
        data_set = self.owner_list_conv()
        for this_type in data_set:
            # Add clusters to objects
            self.out_clusters[this_type] = []
            self.add_clust_to_obj(self.cluster(data_set[this_type]["data"],this_type),
                                  data_set[this_type]["objects"],this_type)


    def add_clust_to_obj(self,cluster_obj,object_list,type):
        """
        Given the output data -> assign an output cluster
        :return:
        """
        for i, cluster in enumerate(cluster_obj.clusters):
            new_cluster = Cluster(cluster_obj.clusters[cluster],type,i)
            self.out_clusters[type].append(new_cluster)
        for i, val in enumerate(cluster_obj.dataClusterId):
            this_clust = [x for x in self.out_clusters[type] if x.cluster == val][0]
            this_obj = object_list[i]
            this_obj.cluster = this_clust
            this_clust.object_list.append(this_obj)

    def write_output(self):
        """
        Write the output of the clustering to JSON or similar
        :return:
        """
        for type in self.out_clusters:
            self.out_clusters[type]

    def owner_list_conv(self):
        """
        Convert an owner list to a Dict of Type - cluster list.
        :return:
        """
        out_d = {}
        for owner in self.owner_list:
            for obj in owner.object_list:
                if obj.object_desc in out_d:
                    out_d[obj.object_desc]["data"].append(obj.value_array)
                    out_d[obj.object_desc]["objects"].append(obj)
                else:
                    out_d[obj.object_desc] = {"data": [obj.value_array],"pks": [obj.uuid]}
        return out_d


class Object(object):

    def __init__(self, value_array, object_desc):
        """
        An Object will have a type and a value arrays
        :param value_array: e.g. x,y,z coords - what is clustered
        :param object_desc: e.g. "Water" or "H-bond acceptor"
        """
        self.value_array = value_array
        self.object_desc = object_desc
        self.uuid = str(uuid.uuid4())
        self.cluster = -1


class Owner(object):

    def __init__(self, object_list, title):
        """
        An Owner will own multiple objects. E.g. PDB 4CUP owns waters.
        :param object_list: the list of objects it owns
        :param title: the title of the object
        """
        self.object_list = object_list
        self.title = title

class Cluster(Object):

    def __init__(self, value_array, object_desc, index):
        """
        An Object will have a type and a value arrays
        :param value_array: e.g. x,y,z coords - what is clustered
        :param object_desc: e.g. "Water" or "H-bond acceptor"
        """
        self.value_array = value_array
        self.object_desc = object_desc
        self.uuid = str(uuid.uuid4())
        self.object_list = []
        self.cluster_ind = index

    def get_size(self):
        return len(self.object_list)