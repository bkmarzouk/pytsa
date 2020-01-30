class bispectrumCfgTemplate:
    
    def __init__(self):
        self.configurations = None
        
    def addConfiguration(self, name, latex, alpha, beta):
        
        assert type(name)  is str, "configuration name must be string"
        assert type(latex) is str, "latex definition must be string"
        assert type(alpha) is float, "alpha parameter must be float"
        assert type(beta)  is float, "alpha parameter must be float"
        
        if self.configurations is None: self.configurations = []
        
        self.configurations.append(
            {
                "name" : name,
                "latex": latex,
                "alpha": alpha,
                "beta" : beta
            }
        )
        
class fieldspaceCfgTemplate:
    
    def __init__(self):
        self.goodExit = None
        self.badExit = None
        
    def addEndRegion(self, fieldNumber, minFieldValue, maxFieldValue):
        
        assert type(fieldNumber)   is int, "field number must be an integer corresponding to the installation"
        assert type(minFieldValue) in [int, float], "minFieldValue must be an integer or float"
        assert type(maxFieldValue) in [int, float], "maxFieldValue must be an integer or float"
        
        if self.goodExit is None: self.goodExit = []
        
        self.goodExit.append(
            {
                "fieldNumber"  : fieldNumber,
                "minFieldValue": minFieldValue,
                "maxFieldValue": maxFieldValue
            }
        )
        
    def addViolatedRegion(self, fieldNumber, minFieldValue, maxFieldValue):
    
        assert type(fieldNumber)   is int, "field number must be an integer corresponding to the installation"
        assert type(minFieldValue) in [int, float], "minFieldValue must be an integer or float"
        assert type(maxFieldValue) in [int, float], "maxFieldValue must be an integer or float"
    
        if self.badExit is None: self.badExit = []
    
        self.badExit.append(
            {
                "fieldNumber"  : fieldNumber,
                "minFieldValue": minFieldValue,
                "maxFieldValue": maxFieldValue
            }
        )
        
class