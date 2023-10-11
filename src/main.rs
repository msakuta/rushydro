mod app;

use app::RusHydroApp;

#[cfg(not(target_arch = "wasm32"))]
fn main() {
    let mut native_options = eframe::NativeOptions::default();

    // We insist to use light theme, because the canvas color is designed to work with light background.
    native_options.follow_system_theme = false;
    native_options.default_theme = eframe::Theme::Light;

    eframe::run_native(
        "RusHydro GUI",
        native_options,
        Box::new(|_cc| Box::new(RusHydroApp::new())),
    )
    .unwrap();
}
