#pragma once

#include <nlohmann/json.hpp>

struct Scales
{
  double length_scale, time_scale, energy_scale;

  void
  from_json(const nlohmann::json &j)
  {
    length_scale = j.at("dimensions").at("length_scale").get<double>();
    time_scale   = j.at("dimensions").at("time_scale").get<double>();
    energy_scale = j.at("dimensions").at("energy_density_scale").get<double>();
  }
};
