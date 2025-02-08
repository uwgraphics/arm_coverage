URDF_PATHS = {  'panda': 'panda.urdf', 
                'ur5': 'ur5.urdf',
                'ur5sander': 'ur5sander.urdf'}
BASE_LINKS = {  'panda': 'world', 
                'ur5': 'base_link',
                'ur5sander': 'base_link'}
TIP_LINKS = {   'panda': 'panda_hand', 
                'ur5': 'tool_tip',
                'ur5sander': 'tool_tip'}
SETTING_FILE_PATHS = {
    'panda': 'example_settings/panda.yaml',
    'ur5': 'example_settings/ur5.yaml',
    'ur5sander': 'example_settings/ur5sander.yaml'
}

LARGE_WEIGHT = 536870911
SEMI_LARGE_WEIGHT = int(1e7-1)
