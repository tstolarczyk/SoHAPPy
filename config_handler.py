import yaml
class ConfigHandler(object):
    
    def __init__(self, config_file):
        """
        Class used to get parameters from
        a configuration file using yaml syntax.

        Parameters
        ----------
        config_file : configuration file (yaml format)
        """
        with open(config_file,'r') as stream:
            self.cfg = yaml.load(stream)

    def get_value(self, param1, param2=None, param3=None):
        """
        Function to get value of variable

        Parameters
        ----------
        param1: section's name
        param2: variable's name or subsection name
        param3: variable's name

        Returns
        -------
        Value parameter (string, float...)
        """
        if param3==None and param2==None:
            return self.cfg[param1]
        elif param3==None:
            return self.cfg[param1][param2]
        else:
            return self.cfg[param1][param2][param3]
