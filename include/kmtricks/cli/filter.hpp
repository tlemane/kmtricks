#pragma once

#include <kmtricks/cli/cli_common.hpp>
#include <kmtricks/cmd/filter.hpp>
#include <kmtricks/config.hpp>

namespace km {

km_options_t filter_cli(std::shared_ptr<bc::Parser<1>> cli, filter_options_t options);

}
