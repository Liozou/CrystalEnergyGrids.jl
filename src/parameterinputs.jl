struct ShootingStarMinimizer
    lb::LoadBalancer{Nothing}
end
function ShootingStarMinimizer()
    lb = LoadBalancer{Nothing}(7) do newmc
        let newmc=newmc
            ignore_args(newmc)
        end
    end
    ShootingStarMinimizer(lb)
end

function (star::ShootingStarMinimizer)(_, _, _, _)
    put!(star.lb, nothing)
    nothing
end

Base.fetch(x::ShootingStarMinimizer) = wait(x.lb)

