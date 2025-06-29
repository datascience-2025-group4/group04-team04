def show_feature_spaces(feature_spaces: dict, regions: list):
    for region in regions:
        fs_show = feature_spaces[region]
        X = fs_show[