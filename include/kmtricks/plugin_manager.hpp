#pragma once
#include <string>
#include <dlfcn.h>
#include <filesystem>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <kmtricks/exceptions.hpp>

namespace fs = std::filesystem;

namespace km {

template<typename P>
class PluginManager
{
private:
  PluginManager() {}

public:
  void init(const std::string& shared_lib_path, const std::string& config, size_t max_size)
  {
    m_max_size = max_size;
    m_config = config;
    m_lib_path = shared_lib_path;
    if (!fs::exists(m_lib_path))
      throw FileNotFoundError(fmt::format("{} not found!", m_lib_path));
    spdlog::info("Load plugin ...");
    load();
  }

  bool use_plugin() const
  {
    return m_enable;
  }

private:
  void load()
  {
    m_handle = dlopen(m_lib_path.c_str(), RTLD_LAZY);

    if (!m_handle)
    {
      throw PluginError(fmt::format("Unable to load shared lib. dlerror: {}", dlerror()));
    }

    int (*use_template)();
    use_template = reinterpret_cast<int(*)()>(dlsym(m_handle, "use_template"));
    const char* dlsym_error = dlerror();
    if (dlsym_error)
      throw PluginError(fmt::format("Unable to load symbol. dlsym_error: {}", dlsym_error));

    bool use_T = use_template();

    load_create(use_T ? m_max_size : 0);

    load_destroy();

    load_name();

    spdlog::info("Plugin '{}' loaded.", m_plugin_name);
    m_enable = true;
  }

  void load_create(size_t max_size)
  {
    m_load_plugin = reinterpret_cast<P*(*)()>(dlsym(m_handle, fmt::format("create{}", max_size).c_str()));

    const char* dlsym_error = dlerror();
    if (dlsym_error)
      throw PluginError(fmt::format("Unable to load symbol. dlsym_error: {}", dlsym_error));
  }

  void load_destroy()
  {
    m_destroy_plugin = reinterpret_cast<void(*)(P*)>(dlsym(m_handle, "destroy"));
    const char* dlsym_error = dlerror();
    if (dlsym_error)
      throw PluginError(fmt::format("Unable to load symbol. dlsym_error: {}", dlsym_error));
  }

  void load_name()
  {
    std::string(*plugin_name)();
    plugin_name = reinterpret_cast<std::string(*)()>(dlsym(m_handle, "plugin_name"));
    const char* dlsym_error = dlerror();
    if (dlsym_error)
      throw PluginError(fmt::format("Unable to load symbol. dlsym_error: {}", dlsym_error));
    m_plugin_name = plugin_name();
  }

public:
  void close()
  {
    if (m_handle)
    {
      dlclose(m_handle);
      m_handle = nullptr;
    }
  }

  ~PluginManager() { close(); }

  static PluginManager<P>& get() { static PluginManager<P> pm; return pm;}

  P* get_plugin()
  {
    P* p = m_load_plugin();
    p->configure(m_config);
    return p;
  }

  void destroy_plugin(P* p)
  {
    m_destroy_plugin(p);
  }

private:
    bool m_enable {false};
    size_t m_max_size;
    std::string m_config;
    std::string m_lib_path;
    bool m_ready {false};
    void* m_handle {nullptr};
    P* (*m_load_plugin)() {nullptr};
    void (*m_destroy_plugin)(P*) {nullptr};
    std::string m_plugin_name;
};

}
