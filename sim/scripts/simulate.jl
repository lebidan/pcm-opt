using BPsim

function main()
    yaml_path, overrides = parse_commandline(ARGS)
    cfg = load_config(yaml_path, overrides)
    run_simulation(cfg)
end

main()
