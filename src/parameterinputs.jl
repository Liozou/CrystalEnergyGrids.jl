struct ShootingStarMinimizer
    lb::LoadBalancer{SimulationStep}
end
function ShootingStarMinimizer()
    lb = LoadBalancer{SimulationStep}(output_restart, 7)
    ShootingStarMinimizer(lb)
end

Base.fetch(x::ShootingStarMinimizer) = wait(x.lb)
